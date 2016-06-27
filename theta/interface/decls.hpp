#ifndef THETA_DECLS_HPP
#define THETA_DECLS_HPP

// this header is included in almost any other header which require the declaration
// of classes. Whenever possible, only this header should be included.

namespace libconfig{
    class Setting;
}

namespace theta{
    class Matrix;
    class DoubleVector;
    class Histogram1D;
    class Histogram1DWithUncertainties;
    class HistogramFunction;
    class Minimizer;
    class Producer;
    class ParameterDependentProducer;
    class ProductsSource;
    class Distribution;
    class Ranges;
    class Database;
    class DatabaseInput;
    class Table;
    class Column;
    class ProductsTable;
    class LogTable;
    class RndInfoTable;
    
    class MCMCStrategy;

    class Configuration;
    class Random;
    class VarIdManager;
    template<typename tag> class VarId;
    struct ParIdTag;
    struct ObsIdTag;
    typedef VarId<ParIdTag> ParId;
    typedef VarId<ObsIdTag> ObsId;
    template<class id_type> class VarIds;
    typedef VarIds<ObsId> ObsIds;
    typedef VarIds<ParId> ParIds;
    class ParValues;
    class PropertyMap;
    
    //class Data;
    template<typename T> class DataT;
    typedef DataT<Histogram1D> Data;
    typedef DataT<Histogram1DWithUncertainties> DataWithUncertainties;
    class Function;
    class Model;
    class DataSource;
    class NLLikelihood;
    
    class Setting;
    class atomic_int;

    class ToyMaker;
    class BufferingProductsSink;
}

#endif
