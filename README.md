# Toeplitz_Matrix_Multiplication
Comparing with time complexity and error in the multiplication of Toeplitz matrix with a general matrix using different algos
## Description
The Toeplitz_multiplication.m scripts finds the inverse of a Toeplitz matrix by using the fact that time complexity of multiplying a circulant matrix with a vector is O(n log n).Toeplitz matrix is converted to circulant matrix by means of two different algorithms. These two algorithms for matrix multiplication are evaluated against Random sampling algorithm for matrix multiplicaiton in terms of time complexity and Error.

## Requirements 
1. MATLAB

## Usage
1. Clone or download the repository to your local machine
2. Open MATLAB and navigate to the directory where you saved the repository.
3. Open the 'Toeplitz_multiplication.m' script in MATLAB.
4. Run the script and observe the output results and plots.

## Results
The script produces plots showing the relative errors in the matrix multiplication found by using three different algorithms and time complexity taken by the three algorithms agianst the time complexity of general matrix multiplication.The plots are drawn for size of the matrix vs errors and time. As matrices are generated randomly for every size of the matrix mean of the error for 100 iterations are considered.

## Applications

Finding the approximate inverse of a Toeplitz matrix has several applications in various fields. Some potential use cases include:

- Signal Processing
- Time Series Analysis
- Numerical Methods
- Linear Systems
- Communication Systems

By providing an optimized and fast algorithm for computing the matrix multiplication involving a Toeplitz matrix, this code offers a valuable tool for researchers, engineers, and practitioners working in these and related fields.

## Contributing
Contributions to this repository are welcome. If you find any issues or have improvements to suggest, feel free to open an issue or submit a pull request.
