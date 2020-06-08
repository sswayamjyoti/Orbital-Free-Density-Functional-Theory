#ifndef SHAPEFUNCTIONS_HH
#define SHAPEFUNCTIONS_HH

#include <dune/common/fvector.hh>

// LinearShapeFunction:
// represents a shape function and provides methods to evaluate the function
// and its gradient
template<class ctype, class rtype, int dim>
class LinearShapeFunction
{
public:
    enum { dimension = dim };

    LinearShapeFunction() : coeff0(0.0), coeff1(0.0) {}

    LinearShapeFunction(rtype coeff0_, const Dune::FieldVector<rtype,dim>& coeff1_)
        : coeff0(coeff0_), coeff1(coeff1_) {}

    void setCoeff(rtype coeff0_, const Dune::FieldVector<rtype,dim>& coeff1_)
    {
        coeff0 = coeff0_;
        coeff1 = coeff1_;
    }

    rtype evaluateFunction(const Dune::FieldVector<ctype,dim>& local) const
    {
        ctype result = coeff0;
        for (int i = 0; i < dim; ++i)
            result += coeff1[i] * local[i];
        return result;
    }

    Dune::FieldVector<rtype,dim>
    evaluateGradient(const Dune::FieldVector<ctype,dim>& local) const
    {
        return coeff1;
    }

private:
    rtype coeff0;
    Dune::FieldVector<rtype,dim> coeff1;
};

// P1ShapeFunctionSet
// initializes one and only one set of LinearShapeFunction
template<class ctype, class rtype, int dim>
class P1ShapeFunctionSet
{
public:
    enum { n = dim + 1 };

    typedef LinearShapeFunction<ctype,rtype,dim> ShapeFunction;
    typedef rtype resulttype;

    // get the only instance of this class
    static const P1ShapeFunctionSet& instance()
    {
        static const P1ShapeFunctionSet sfs;
        return sfs;
    }

    const ShapeFunction& operator[](int i) const
    {
        if (!i)
            return f0;
        else
            return f1[i - 1];
    }

private:
    // private constructor prevents additional instances
    P1ShapeFunctionSet()
    {
        Dune::FieldVector<rtype,dim> e(-1.0);
        f0.setCoeff(1.0, e);
        for (int i = 0; i < dim; ++i)
        {
            Dune::FieldVector<rtype,dim> e(0.0);
            e[i] = 1.0;
            f1[i].setCoeff(0.0, e);
        }
    }

    P1ShapeFunctionSet(const P1ShapeFunctionSet& other)
    {}

    ShapeFunction f0;
    ShapeFunction f1[dim];
};

#endif
