
// Function for sequence alignment with affine gap penalty
int sequenceAlignment(const std::string& x, const std::string& y, int gapOpen, int gapExtend) {
    int m = x.length();
    int n = y.length();

    // Initialize matrices for dynamic programming
    std::vector<std::vector<int>> M(m + 1, std::vector<int>(n + 1, 0));  // Match/mismatch matrix
    std::vector<std::vector<int>> A(2, std::vector<int>(n + 1, 0));  // Row A
    std::vector<std::vector<int>> B(2, std::vector<int>(n + 1, 0));  // Row B
    std::vector<std::vector<int>> C(m + 1, std::vector<int>(n + 1, 0));  // Substitution matrix

    // Initialization
    for (int i = 1; i <= m; ++i) {
        M[i][0] = -gapOpen - (i - 1) * gapExtend;
    }

    for (int j = 1; j <= n; ++j) {
        M[0][j] = -gapOpen - (j - 1) * gapExtend;
    }

    // Dynamic programming
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            // Compute substitution matrix C
            C[i][j] = (x[i - 1] == y[j - 1]) ? 0 : 1;

            // Update row A with gap extension
            std::transform(M[i - 1].begin(), M[i - 1].end(), A[0].begin(), [gapExtend](int val) {
                return val - gapExtend;
            });

            // Update row B with gap opening
            copy(A[0].begin(), A[0].end(), B[0].begin());
            transform(B[0].begin() + 1, B[0].end(), B[1].begin(), [gapOpen](int val) {
                return val - gapOpen;
            });

            // Update row A with substitution matrix
            transform(M[i - 1].begin() + 1, M[i - 1].end(), C[i].begin() + 1, A[1].begin(), std::plus<int>());

            // Element-wise minimum of A, B, and C
            transform(A[1].begin(), A[1].end(), B[1].begin(), M[i].begin() + 1, [C](int a, int b) {
                return minimum(a, b, C[i][j]);
            });
        }
    }

    // The final result is stored in M[m][n]
    return M[m][n];
}