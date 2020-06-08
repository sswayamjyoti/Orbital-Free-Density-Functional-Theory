inline double POT_EL_ATOMS(const std::list<ATOM> atoms, const Dune::FieldVector<double,3>& pos_qp)
{
       double result = 0.0;
       Dune::FieldVector<double,3> dist;
       for (auto ait= atoms.begin(); ait != atoms.end(); ++ait)
       {
         dist = ait->pos;
         dist -= pos_qp;
         result += (ait->charge / dist.two_norm());
       }
       return result;
}; 
