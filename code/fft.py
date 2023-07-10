import math
j = complex(0, 1)

def round_complex_number(c):
    return round(c.real, 10) + round(c.imag, 10) * 1j

def round_complex_array(C):
    return [round_complex_number(c) for c in C]

# Naive Fourier Transform
def naive_discrete_fourier_transform(x):
    N = len(x)
    X = [None] * N
    # X_k = sum from n=0..N-1 of x_n * e^(-j*2*pi/N*k*n)
    for k in range(N):
        X[k] = sum(x[n] * math.e ** (-j * 2 * math.pi / N * k * n) for n in range(N))
    return round_complex_array(X)

def inverse_naive_discrete_fourier_transform(X):
    N = len(X)
    x = [None] * N
    # x_n = 1/N * sum from k=0..N-1 of X_k * e^(j*2*pi/N*k*n)
    for n in range(N):
        x[n] = 1 / N * sum(X[k] * math.e ** (j * 2 * math.pi / N * k * n) for k in range(N))
    return round_complex_array(x)

x = [0, 1, 2, 3]
X = naive_discrete_fourier_transform(x)
print(X)
x = inverse_naive_discrete_fourier_transform(X)
print(x)

# Fast Fourier Transform
def fast_fourier_transform_recursive(x):
    N = len(x)
    if N == 1:
        return x
    # principle Nth root of unity (negative version)
    w_N = math.e ** (-j * 2 * math.pi / N)
    # current twiddle factor
    w = 1
    # compute the Fourier Transform of the even and odd indices
    x_even = [x[n] for n in range(0, N, 2)]
    x_odd = [x[n] for n in range(1, N, 2)]
    X_even = fast_fourier_transform_recursive(x_even)
    X_odd = fast_fourier_transform_recursive(x_odd)
    X = [None] * N
    for k in range(N // 2):
        # X[k] and X[k + N // 2] are computed using the same even and odd coefficients
        # by just reversing the sign of the twiddle factor of the odd component
        X[k] = X_even[k] + w * X_odd[k]
        X[k + N // 2] = X_even[k] - w * X_odd[k]
        w = w * w_N
    return X

def fast_fourier_transform(x):
    N = len(x)
    assert N >= 1 and (N & (N-1)) == 0
    X = fast_fourier_transform_recursive(x)
    return round_complex_array(X)

def inverse_fast_fourier_transform_recursive(X):
    N = len(X)
    if N == 1:
        return X
    # principle Nth root of unity (positive version)
    w_N = math.e ** (j * 2 * math.pi / N)
    # current twiddle factor
    w = 1
    X_even = [X[k] for k in range(0, N, 2)]
    X_odd = [X[k] for k in range(1, N, 2)]
    x_even = inverse_fast_fourier_transform_recursive(X_even)
    x_odd = inverse_fast_fourier_transform_recursive(X_odd)
    x = [None] * N
    for n in range(N // 2):
        # x[n] and x[n + N // 2] are computed using the same even and odd coefficients
        # by just reversing the sign of the twiddle factor of the odd component
        x[n] = x_even[n] + w * x_odd[n]
        x[n + N // 2] = x_even[n] - w * x_odd[n]
        w = w * w_N
    return x

def inverse_fast_fourier_transform(X):
    N = len(X)
    assert N >= 1 and (N & (N-1)) == 0
    x = inverse_fast_fourier_transform_recursive(X)
    # normalize by dividing by N
    x_normalized = [x_n / N for x_n in x]
    return round_complex_array(x_normalized)

x = [0, 1, 2, 3]
X = fast_fourier_transform(x)
print(X)
x = inverse_fast_fourier_transform(X)
print(x)

# Polynomial Multiplication
def polynomial_multiplication(A, B):
    a_length, b_length = len(A), len(B)
    N = max(a_length, b_length)
    # round N to next largest power of two
    N = 2 ** math.ceil(math.log(N, 2))
    # zerofill coefficient arrays to 2N
    A = A + [0] * (2 * N - len(A))
    B = B + [0] * (2 * N - len(B))
    # use FFT to get point-value representation of A and B
    A_pv = fast_fourier_transform(A)
    B_pv = fast_fourier_transform(B)
    # multiply point-values to get point-value representation of C
    C_pv = [a_pv * b_pv for a_pv, b_pv in zip(A_pv, B_pv)]
    # use IFFT to get coefficient representation of C
    C = inverse_fast_fourier_transform(C_pv)
    # round and concatenate to correct length
    return [round(c.real) for c in C[:a_length + b_length]]

A = [1, 2, 3, 4]
B = [5, 6, 7, 8]
print(polynomial_multiplication(A, B))

"""
def solve():
    T = int(input())
    for t in range(T):
        input()
        A = list(map(int, input().split()))
        B = list(map(int, input().split()))
        print(" ".join(map(str, polynomial_multiplication(A, B))))
solve()
"""
