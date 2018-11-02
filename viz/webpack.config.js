var path = require('path');

module.exports = {
  entry: ['@babel/polyfill', './js/main'],
  output: {
    filename: '../dist.js'
  },
  devtool: 'source-map',
  module: {
    rules: [{
      test: /\.js$/,
      exclude: /node_modules/,
      use: {
        loader: 'babel-loader',
        options: {
          presets: ['@babel/preset-env', '@babel/preset-react'],
          plugins: ['@babel/plugin-proposal-class-properties']
        }
      }
    }]
  },
  resolve: {
    extensions: ['.js']
  }
};