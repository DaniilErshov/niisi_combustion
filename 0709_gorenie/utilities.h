#pragma once


#include <numeric>

    //! @addtogroup mathTemplates
    //! @{

    //! Templated Inner product of two vectors of length 4.
    /*!
     * If either @e x or @e y has length greater than 4, only the first 4 elements
     * will be used.
     *
     * @param x   first reference to the templated class V
     * @param y   second reference to the templated class V
     * @return This class returns a hard-coded type, double.
     */
    template<class V>
    inline double dot4(const V& x, const V& y)
    {
        return x[0] * y[0] + x[1] * y[1] + x[2] * y[2] + x[3] * y[3];
    }

    //! Templated Inner product of two vectors of length 5
    /*!
     * If either @e x or @e y has length greater than 4, only the first 4 elements
     * will be used.
     *
     * @param x   first reference to the templated class V
     * @param y   second reference to the templated class V
     * @return This class returns a hard-coded type, double.
     */
    template<class V>
    inline double dot5(const V& x, const V& y)
    {
        return x[0] * y[0] + x[1] * y[1] + x[2] * y[2] + x[3] * y[3] +
            x[4] * y[4];
    }

    //! Function that calculates a templated inner product.
    /*!
     * This inner product is templated twice. The output variable is hard coded
     * to return a double.
     *
     * template<class InputIter, class InputIter2>
     *
     * @code
     *     double x[8], y[8];
     *     double dsum = dot<double *,double *>(x, &x+7, y);
     * @endcode
     *
     * @param x_begin  Iterator pointing to the beginning, belonging to the
     *                 iterator class InputIter.
     * @param x_end    Iterator pointing to the end, belonging to the
     *                 iterator class InputIter.
     * @param y_begin Iterator pointing to the beginning of y, belonging to the
     *               iterator class InputIter2.
     * @return The return is hard-coded to return a double.
     */
    template<class InputIter, class InputIter2>
    inline double dot(InputIter x_begin, InputIter x_end, InputIter2 y_begin)
    {
        return std::inner_product(x_begin, x_end, y_begin, 0.0);
    }

    //! Multiply elements of an array by a scale factor.
    /*!
     * @code
     * vector<double> in(8, 1.0), out(8);
     * scale(in.begin(), in.end(), out.begin(), factor);
     * @endcode
     *
     * @param begin  Iterator pointing to the beginning, belonging to the
     *               iterator class InputIter.
     * @param end    Iterator pointing to the end, belonging to the
     *               iterator class InputIter.
     * @param out    Iterator pointing to the beginning of out, belonging to the
     *               iterator class OutputIter. This is the output variable
     *               for this routine.
     * @param scale_factor  input scale factor belonging to the class S.
     */
    //! Templated evaluation of a polynomial of order 6
    /*!
     * @param x   Value of the independent variable - First template parameter
     * @param c   Pointer to the polynomial - Second template parameter
     */
    template<class D, class R>
    R poly6(D x, R* c)
    {
        return ((((((c[6] * x + c[5]) * x + c[4]) * x + c[3]) * x +
            c[2]) * x + c[1]) * x + c[0]);
    }

    //! Templated evaluation of a polynomial of order 8
    /*!
     * @param x   Value of the independent variable - First template parameter
     * @param c   Pointer to the polynomial - Second template parameter
     */
    template<class D, class R>
    R poly8(D x, R* c)
    {
        return ((((((((c[8] * x + c[7]) * x + c[6]) * x + c[5]) * x + c[4]) * x + c[3]) * x +
            c[2]) * x + c[1]) * x + c[0]);
    }

    //! Templated evaluation of a polynomial of order 5
    /*!
     * @param x   Value of the independent variable - First template parameter
     * @param c   Pointer to the polynomial - Second template parameter
     */
    template<class D, class R>
    R poly5(D x, R* c)
    {
        return (((((c[5] * x + c[4]) * x + c[3]) * x +
            c[2]) * x + c[1]) * x + c[0]);
    }

    //! Evaluates a polynomial of order 4.
    /*!
     * @param x   Value of the independent variable.
     * @param c   Pointer to the polynomial coefficient array.
     */
    template<class D, class R>
    R poly4(D x, R* c)
    {
        return ((((c[4] * x + c[3]) * x +
            c[2]) * x + c[1]) * x + c[0]);
    }

    //! Templated evaluation of a polynomial of order 3
    /*!
     * @param x   Value of the independent variable - First template parameter
     * @param c   Pointer to the polynomial - Second template parameter
     */
    template<class D, class R>
    R poly3(D x, R* c)
    {
        return (((c[3] * x + c[2]) * x + c[1]) * x + c[0]);
    }