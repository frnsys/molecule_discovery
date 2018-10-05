# system installation
DIR=$(pwd)

cd /tmp
sudo apt install unzip
wget https://github.com/bazelbuild/bazel/releases/download/0.16.0/bazel-0.16.0-installer-linux-x86_64.sh
chmod +x bazel-0.16.0-installer-linux-x86_64.sh
sudo ./bazel-0.16.0-installer-linux-x86_64.sh
sudo pip3 install -U keras

git clone https://github.com/tensorflow/tensorflow
cd tensorflow
./configure

# optionally add --config=cuda for CUDA
bazel build -c opt --copt=-mavx512f --copt=-mavx --copt=-mavx2 --copt=-mfma --copt=-mfpmath=both -k //tensorflow/tools/pip_package:build_pip_package
bazel-bin/tensorflow/tools/pip_package/build_pip_package /tmp/tensorflow_pkg
sudo pip3 install /tmp/tensorflow_pkg/tensorflow-*-linux_x86_64.whl

cd /tmp
git clone https://github.com/abseil/abseil-cpp.git
sudo mv abseil-cpp/absl /usr/local/lib/python3.6/dist-packages/tensorflow/include/absl

cd $DIR
git clone https://github.com/tensorflow/models.git
cd models/tutorials/embedding
TF_CFLAGS=( $(python3 -c 'import tensorflow as tf; print(" ".join(tf.sysconfig.get_compile_flags()))') )
TF_LFLAGS=( $(python3 -c 'import tensorflow as tf; print(" ".join(tf.sysconfig.get_link_flags()))') )
g++ -std=c++11 -shared word2vec_ops.cc word2vec_kernels.cc -o word2vec_ops.so -fPIC ${TF_CFLAGS[@]} ${TF_LFLAGS[@]} -O2
cp "$DIR/word2vec.py" word2vec.py
python3 word2vec.py --save_path="$DIR/w2v" --train_data="$DIR/w2v.txt"
