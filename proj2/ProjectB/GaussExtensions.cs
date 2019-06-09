using System;
using Core;

namespace ProjectB
{
    public static class GaussExtensions
    {
        /// <summary>
        /// This function computes the elementary row replacement operation on
        /// the given matrix.
        /// </summary>
        ///
        /// <remarks>
        /// Note that we add the row (as in the lectures) instead of subtracting
        /// the row (as in the textbook).
        /// </remarks>
        ///
        /// <param name="a">
        /// An N-by-M matrix to perform the elementary row operation on.
        /// </param>
        /// <param name="i">
        /// The index of the row to replace.
        /// </param>
        /// <param name="m">
        /// The multiplum of row j to add to row i.
        /// </param>
        /// <param name="j">
        /// The index of the row to replace with.
        /// </param>
        ///
        /// <returns>
        /// The resulting N-by-M matrix after having performed the elementary
        /// row operation.
        /// </returns>
        public static Matrix ElementaryRowReplacement(
            this Matrix a, int i, double m, int j) {
            for (int n = 0; n < a.N_Cols; n++) {
                a[i, n] += a[j, n] * m;
            }
            return a;
        }

        /// <summary>
        /// This function computes the elementary row interchange operation on
        /// the given matrix.
        /// </summary>
        ///
        /// <param name="a">
        /// An N-by-M matrix to perform the elementary row operation on.
        /// </param>
        /// <param name="i">
        /// The index of the first row of the rows to interchange.
        /// </param>
        /// <param name="j">
        /// The index of the second row of the rows to interchange.
        /// </param>
        ///
        /// <returns>
        /// The resulting N-by-M matrix after having performed the elementary
        /// row operation.
        /// </returns>
        public static Matrix ElementaryRowInterchange(
            this Matrix a, int i, int j)
        {
            for(var n =0; n < a.N_Cols; n++) {
                var tmp = a[i, n];
                a[i, n] = a[j, n];
                a[j, n] = tmp;

            }

            return a;
        }

        /// <summary>
        /// This function computes the elementary row scaling operation on the
        /// given matrix.
        /// </summary>
        ///
        /// <param name="a">
        /// An N-by-M matrix to perform the elementary row operation on.
        /// </param>
        /// <param name="i">The index of the row to scale.</param>
        /// <param name="c">The value to scale the row by.</param>
        ///
        /// <returns>
        /// The resulting N-by-M matrix after having performed the elementary
        /// row operation.
        /// </returns>
        public static Matrix ElementaryRowScaling(
            this Matrix a, int i, double c)
        {
            for (int n = 0; n < a.N_Cols; n++) {
                a[i, n] *= c;
            }

            return a;
        }

        static bool is_non_zero_col(Matrix a, int j) {
            for (var m = 0; m <= a.M_Rows; m++) {
                if (a[j, m] == 1)
                    return true;
            }

            return false;
        }
        
        static bool is_zero_row(Matrix a, int row){
            for (int i = 0; i < a.N_Cols; i++) {
                if (!is_zero_with_margin(a[row, i], 1e-8)) {
                    return true;
                }
            }

            return false;
        }

        static bool is_zero_with_margin(double f, double margin) {
            return Math.Abs(f) < margin;
        }
        
        static int find_pivot_in_col(Vector v, int start, double margin) {
            for (int i = start; i < v.Size; i++) {
                // If the number is zero
                if (!is_zero_with_margin(v[i], margin)) {
                    return i;
                }
            }
            return -1;
        }

        /// <summary>
        /// This function executes the forward reduction algorithm provided in
        /// the assignment text to achieve row Echelon form of a given
        /// augmented matrix.
        /// </summary>
        ///
        /// <param name="a">
        /// An N-by-M augmented matrix.
        /// </param>
        ///
        /// <returns>
        /// An N-by-M matrix that is the row Echelon form.
        /// </returns>
        public static Matrix ForwardReduction(this Matrix a)
        {
            // We suggest this tolerance number for floating point comparisons.
            var tol = 1e-8;

            var p = 0;
            for (var j = 0; j < a.N_Cols; j++) {
                // Find(and maybe swap) pivot element
                var i = find_pivot_in_col(a.Column(j), p, tol);
                if (i == -1) { // Not a pivot col
                    continue;
                }

                if (i != p) { // index is in correct row, swap
                    a.ElementaryRowInterchange(i, p);
                }

                var pivot = a[p, j]; // Now scale the row to 1
                if (pivot <= 1 - tol || pivot >= 1 + tol) {
                    a.ElementaryRowScaling(p, Math.Pow(pivot, -1));
                }

                // And now make 0 all the way down
                for (var k = p + 1; k < a.M_Rows; k++) {
                    if (!is_zero_with_margin(a[k, j], tol)) {
                        a.ElementaryRowReplacement(k, -a[k, j], p);
                    }
                }

                a.ElementaryRowScaling(p, pivot);

                p++;
                if (p == a.M_Rows) {
                    break;
                }
            }
            
            return a;
        }

        /// <summary>
        /// This function executes the backward reduction algorithm provided in
        /// the assignment text given an augmented matrix in row Echelon form.
        /// </summary>
        ///
        /// <param name="a">
        /// An N-by-M augmented matrix in row Echelon form.
        /// </param>
        ///
        /// <returns>
        /// The resulting N-by-M matrix after executing the algorithm.
        /// </returns>
        public static Matrix BackwardReduction(this Matrix a)
        {
            // We suggest this tolerance number for floating point comparisons.
            var tol = 1e-8;

            for (var p = a.M_Rows - 1; p >= 0; p--) {
                // Find pivot in row
                var j = 0;
                while(j < a.N_Cols) {
                    if (!is_zero_with_margin(a[p, j], tol)) {
                        break;
                    }
                    j++;
                }

                if (j == a.N_Cols) {// Zero row
                    continue;
                }
                // Scale row to 1
                a.ElementaryRowScaling(p, Math.Pow(a[p, j], -1));
                for (var i = p - 1; i >= 0; i--) {
                    a.ElementaryRowReplacement(i, -a[i, j], p);
                }
            }

            return a;

        }

        /// <summary>
        /// This function performs Gauss elimination of a linear system
        /// given in matrix form by a coefficient matrix and a right hand side
        /// vector. It is assumed that the corresponding linear system is
        /// consistent and has exactly one solution.
        /// </summary>
        ///
        /// <remarks>
        /// Hint: Combine ForwardReduction and BackwardReduction.
        /// </remarks>
        ///
        /// <param name="a">An N-by-M matrix.</param>
        /// <param name="b">An N-size vector.</param>
        ///
        /// <returns>The M-sized vector x such that a * x = b.</returns>
        public static Vector GaussElimination(this Matrix a, Vector b) {
            return a.AugmentRight(b).ForwardReduction().BackwardReduction().Column(a.N_Cols);
        }
    }

    public static class BasicExtensions
    {
        /// <summary>
        /// This function creates an augmented matrix given a matrix 'a' and a
        /// right-hand side vector 'v'.
        /// </summary>
        ///
        /// <remarks>
        /// See page 12 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <param name="a">An M-by-N matrix.</param>
        /// <param name="v">An M-size vector.</param>
        ///
        /// <returns>The M-by-(N + 1) augmented matrix [a | v].</returns>
        public static Matrix AugmentRight(this Matrix a, Vector v)
        {
            var mRows = a.M_Rows;
            var nCols = a.N_Cols;

            var retval = new double[mRows, nCols + 1]; // 0-initialized

            for (var i = 0; i < mRows; i++)
            {
                for (var j = 0; j < nCols; j++)
                {
                    retval[i, j] = a[i, j];
                }
                retval[i, nCols] = v[i];
            }

            return new Matrix(retval);
        }
    }
}
