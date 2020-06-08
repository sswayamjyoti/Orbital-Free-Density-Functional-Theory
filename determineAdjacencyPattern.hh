
template<class GV> 
void P1Elements<GV>::determineAdjacencyPattern()
{
    const int N = gv.size(dim);
    adjacencyPattern.resize(N);

    const LeafIndexSet& set = gv.indexSet();
    const LeafIterator itend = gv.template end<0>();

    for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
    {
        Dune::GeometryType gt = it->type();
        const Dune::template GenericReferenceElement<ctype,dim> &ref =
            Dune::GenericReferenceElements<ctype,dim>::general(gt);

        // traverse all codim-1-entities of the current element and store all
        // pairs of vertices in adjacencyPattern
        const IntersectionIterator isend = gv.iend(*it);
        for (IntersectionIterator is = gv.ibegin(*it) ; is != isend ; ++is)
        {
            int vertexsize = ref.size(is->indexInInside(),1,dim);
            for (int i=0; i < vertexsize; i++)
            {
                int indexi = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,i,dim),dim);
                for (int j=0; j < vertexsize; j++)
                {
                    int indexj = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,j,dim),dim);
                    adjacencyPattern[indexi].insert(indexj);
                }
            }
        }
    }
} 
