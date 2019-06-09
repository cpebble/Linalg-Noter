using System;
using Core;

namespace ProjectC
{
    public static class AdvancedExtensions
    {
        /// <summary>
        /// This function creates the square submatrix given a square matrix as
        /// well as row and column indices to remove from it.
        /// </summary>
        ///
        /// <remarks>
        /// See page 246-247 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <param name="a">An N-by-N matrix.</param>
        /// <param name="i">The index of the row to remove.</param>
        /// <param name="j">The index of the column to remove.</param>
        ///
        /// <returns>The resulting (N - 1)-by-(N - 1) submatrix.</returns>
        public static Matrix SquareSubMatrix(this Matrix a, int i, int j)
        {
            var ret = new Matrix(a.N_Cols - 1, a.N_Cols - 1);

            for (int m = 0; m < ret.M_Rows; m++) {
                for (int n = 0; n < ret.N_Cols; n++) {
                    int[] target = {m,n};
                    if (m >= i) {
                        target[0]++;
                    }

                    if (n >= j) {
                        target[1]++;
                    }

                    ret[m, n] = a[target[0], target[1]];

                }
            }

            return ret;
        }

        /// <summary>
        /// This function computes the determinant of a given square matrix.
        /// </summary>
        ///
        /// <remarks>
        /// See page 247 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <remarks>
        /// Hint: Use SquareSubMatrix.
        /// </remarks>
        ///
        /// <param name="a">An N-by-N matrix.</param>
        ///
        /// <returns>The determinant of the matrix</returns>
        public static double Determinant(this Matrix a) {
            if (a.N_Cols == 1) { // Base case is 1x1 matrix
                return a[0, 0];
            }
            var ret = 0.0; // Running product of determinant
            var sign = 1; // Flips the sign of determinant
            for (int i = 0; i < a.N_Cols; i++) {
                ret += a[0, i] * Determinant(SquareSubMatrix(a, 0, i)) * sign;
                sign = -sign; // We add with alternate signs
            }

            return ret;
        }

        /// <summary>
        /// Calculate ||v||
        /// </summary>
        /// <param name="v">A vector</param>
        /// <returns>The norm of vector v</returns>
        public static double VectorNorm(this Vector v) {
            var ret = 0.0;
            for (int i = 0; i < v.Size; i++) {
                ret += (Math.Pow(v[i], 2));
            }

            return Math.Sqrt(ret);
        }

        /// <summary>
        /// The sum of a vector
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static double VectorSum(this Vector v) {
            var sum = 0.0;
            for (var i = 0; i < v.Size; i++)
                sum += v[i];
            return sum;
        }

        public static double DotProduct(Vector v, Vector u) {
            var ret = 0.0;
            for (int i = 0; i < v.Size; i++) {
                ret += v[i] * u[i];
            }

            return ret;
        }
        /// <summary>
        /// This function computes the Gram-Schmidt process on a given matrix.
        /// </summary>
        ///
        /// <remarks>
        /// See page 229 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <param name="a">
        /// An M-by-N matrix. All columns are implicitly assumed linear
        /// independent.
        /// </param>
        ///
        /// <returns>
        /// A tuple (Q,R) where Q is a M-by-N orthonormal matrix and R is an
        /// N-by-N upper triangular matrix.
        /// </returns>
        public static Tuple<Matrix, Matrix> GramSchmidt(this Matrix a)
        {
            // We suggest this tolerance number for floating point comparisons.
            var tol = 1e-8;
            var Q = new Matrix(a.M_Rows, a.N_Cols);
            var R = new Matrix(a.N_Cols, a.N_Cols);
            for (int j = 0; j < a.N_Cols; j++) {
                var qj = a.Column(j);
                for (int i = 0; i < j; i++) {
                    R[i, j] = DotProduct(Q.Column(i), a.Column(j));
                    qj -= R[i, j] * Q.Column(i);
                }

                var sum = qj.VectorSum();
                if (sum < tol && sum > -tol) {
                    continue;
                }
                R[j, j] = qj.VectorNorm();
                for (int k = 0;k < Q.M_Rows;k++) {
                    Q[k, j] = qj[k] / R[j, j];
                }

            }
            return new Tuple<Matrix, Matrix>(Q,R);
        }
    }
}
