import os
from atc import ATCModel, load_atc_data

# Train the ATC code prediction model
atc_data, atc_idx = load_atc_data(level=3)
if not os.path.exists('data/atc/checkpoint'):
    print('Creating new')
    model = ATCModel(n_features=2048, n_classes=len(atc_idx), learning_rate=0.001, layers=[128,128,128])
else:
    print('Loading')
    model = ATCModel.load('data/atc')
model.train(atc_data, epochs=4000)
save_path = model.save('data/atc')