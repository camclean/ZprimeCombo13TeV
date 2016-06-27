#include "interface/random.hpp"
#include "interface/variables.hpp"

#include <boost/test/unit_test.hpp>

using namespace theta;
using namespace std;

BOOST_AUTO_TEST_SUITE(variables)

struct varidtest{
    VarIdManager vm;
    const ParId par0;
    varidtest(): par0(vm.create_par_id("par0")){
    }
};


BOOST_AUTO_TEST_CASE(basic){
    varidtest v;
    BOOST_CHECK(v.vm.get_name(v.par0)=="par0");
    BOOST_CHECK(v.vm.get_par_id("par0")==v.par0);
}

BOOST_AUTO_TEST_CASE(get_par_ids){
    varidtest v;
    ParIds ids;
    ids.insert(v.par0);
    BOOST_CHECK(ids==v.vm.get_all_parameters());
    
    ParId par1 = v.vm.create_par_id("par1");
    BOOST_REQUIRE(par1!=v.par0);
    BOOST_CHECK(not (ids==v.vm.get_all_parameters()));
    
    ids.insert(par1);
    BOOST_CHECK(ids==v.vm.get_all_parameters());
}

BOOST_AUTO_TEST_CASE(par_exceptions){
    varidtest v;
    
    //request variables not there:
    bool ex = false;
    try{
        v.vm.get_par_id("var1");
    }
    catch(invalid_argument &){
        ex = true;
    }
    BOOST_REQUIRE(ex);
    
    //create parameter already there
    ex = false;
    try{
        v.vm.create_par_id("par0");
    }
    catch(invalid_argument &){
        ex = true;
    }
    BOOST_REQUIRE(ex);
    
    ParId par1 = v.vm.create_par_id("par1");
    BOOST_REQUIRE(v.par0!=par1);
    ParIds all_ids;
    all_ids.insert(par1);
    all_ids.insert(v.par0);
    BOOST_REQUIRE(all_ids == v.vm.get_all_parameters());
}

BOOST_AUTO_TEST_CASE(basic_obs){
    VarIdManager vm;
    ObsId obs0 = vm.create_obs_id("obs0", 100, -0.2, 0.8);
    BOOST_CHECK(vm.get_obs_id("obs0")==obs0);
    BOOST_CHECK(vm.get_nbins(obs0)==100);
    BOOST_CHECK(vm.get_range(obs0).first==-0.2);
    BOOST_CHECK(vm.get_range(obs0).second==0.8);
}

BOOST_AUTO_TEST_CASE(all_obs){
    VarIdManager vm;
    ObsIds ids;
    BOOST_REQUIRE(vm.get_all_observables()==ids);
    ObsId obs0 = vm.create_obs_id("obs0", 100, -0.8, 1.2);
    BOOST_REQUIRE(not (vm.get_all_observables()==ids));
    ids.insert(obs0);
    BOOST_REQUIRE(vm.get_all_observables()==ids);
}

BOOST_AUTO_TEST_CASE(exceptions_obs){
    VarIdManager vm;
    ObsIds ids;
    bool ex = false;
    try{
        vm.create_obs_id("obs0", 0, -1, 1);
    }
    catch(invalid_argument &){
        ex = true;
    }
    BOOST_CHECK(ex);
    BOOST_REQUIRE(vm.get_all_observables()==ids);
    
    ex = false;
    try{
        vm.create_obs_id("obs0", 100, 1, -1);
    }
    catch(invalid_argument &){
        ex = true;
    }
    BOOST_CHECK(ex);
    BOOST_REQUIRE(vm.get_all_observables()==ids);
}


BOOST_AUTO_TEST_CASE(parvalues_basic){
    ParValues vv;
    VarIdManager vm;
    ParId v0 = vm.create_par_id("v0");
    BOOST_REQUIRE(not vv.contains(v0));
    vv.set(v0, -7.0);
    BOOST_REQUIRE(vv.get(v0)==-7.0);
    BOOST_REQUIRE(vv.contains(v0));
    ParId v1 = vm.create_par_id("v1");
    bool ex = false;
    try{
        vv.get(v1);
    }
    catch(invalid_argument &){
        ex = true;
    }
    BOOST_REQUIRE(ex);
}

//as ParValues uses an internally growing vector, test with many parameters ...
BOOST_AUTO_TEST_CASE(parvalues_many){
    ParValues vv;
    VarIdManager vm;
    vector<ParId> parameters;
    ParValues values;
    for(size_t i=0; i<100; ++i){
        stringstream name;
        name << "parameter" << i;
        parameters.push_back(vm.create_par_id(name.str()));
        values.set(parameters.back(), 0.1*i + i*i);
    }
    for(size_t i=0; i<100; ++i){
        BOOST_REQUIRE_EQUAL(values.get(parameters[i]), 0.1*i + i*i);
    }
    
    ParValues values_sparse;
    values_sparse.set(parameters[90], 99.0);
    for(size_t i=0; i<100; ++i){
        if(i==90){
            BOOST_CHECK(values_sparse.get(parameters[i])==99.0);
        }
        else{
            BOOST_REQUIRE(not values_sparse.contains(parameters[i]));
        }
    }
    
    values.set(values_sparse);
    for(size_t i=0; i<100; ++i){
        if(i!=90){
            BOOST_REQUIRE_EQUAL(values.get(parameters[i]), 0.1*i + i*i);
        }
        else{
            BOOST_CHECK(values.get(parameters[i]) == values_sparse.get(parameters[i]));
        }
    }
}

BOOST_AUTO_TEST_CASE(par_ids){
    VarIdManager vm;
    vector<ParId> vpids;
    for(int i=0; i<200; ++i){
        stringstream ss;
        ss << "p" << i;
        vpids.push_back(vm.create_par_id(ss.str()));
    }
    ParIds pids;
    BOOST_CHECK(pids.size() == 0);
    BOOST_CHECK(pids.begin() == pids.end());

    // 1 low id:
    pids.insert(vpids[0]);
    BOOST_CHECK(pids.size()==1);
    BOOST_CHECK(*pids.begin() == vpids[0]);
    ParIds::const_iterator it = pids.begin();
    ++it;
    BOOST_CHECK(it == pids.end());

    // two low ids:
    pids = ParIds();
    pids.insert(vpids[10]);
    pids.insert(vpids[20]);
    BOOST_CHECK(pids.size()==2);
    it = pids.begin();
    BOOST_CHECK(it!=pids.end());
    BOOST_CHECK(*it == vpids[10]);
    ++it;
    BOOST_CHECK(it!=pids.end());
    BOOST_CHECK(*it == vpids[20]);
    ++it;
    BOOST_CHECK(it == pids.end());
    // an additional large one:
    pids.insert(vpids[101]);
    BOOST_CHECK(pids.size()==3);
    it = pids.begin();
    BOOST_CHECK(it!=pids.end());
    BOOST_CHECK(*it == vpids[10]);
    ++it;
    BOOST_CHECK(it!=pids.end());
    BOOST_CHECK(*it == vpids[20]);
    ++it;
    BOOST_CHECK(it!=pids.end());
    BOOST_CHECK(*it == vpids[101]);
    ++it;
    BOOST_CHECK(it == pids.end());



    // 1 large id:
    pids = ParIds();
    pids.insert(vpids[100]);
    BOOST_CHECK(pids.size()==1);
    //ParId id = *pids.begin();
    //cout << "100: " << id << endl;
    BOOST_CHECK(*pids.begin() == vpids[100]);
    it = pids.begin();
    BOOST_CHECK(it!=pids.end());
    ++it;
    BOOST_CHECK(it == pids.end());
    // double insert:
    pids.insert(vpids[100]);
    BOOST_CHECK(pids.size()==1);


    // border ids:
    pids = ParIds();
    pids.insert(vpids[63]);
    BOOST_CHECK(pids.size()==1);
    it = pids.begin();
    BOOST_CHECK(it!=pids.end());
    BOOST_CHECK(*it == vpids[63]);
    ++it;
    BOOST_CHECK(it==pids.end());

    pids = ParIds();
    pids.insert(vpids[64]);
    BOOST_CHECK(pids.size()==1);
    it = pids.begin();
    BOOST_CHECK(it!=pids.end());
    BOOST_CHECK(*it == vpids[64]);
    ++it;
    BOOST_CHECK(it==pids.end());

    // almost all ids:
    for(int i=0; i<200; ++i){
        if(i!=10 && i!=102 && i!=198)
            pids.insert(vpids[i]);
    }
    BOOST_CHECK(pids.size()==197);
    int i=0;
    for(it=pids.begin(); it!=pids.end(); ++it){
        if(i!=10 && i!=102 && i!=198){
            BOOST_CHECK(*it == vpids[i]);
            ++i;
        }
    }
    BOOST_CHECK(it==pids.end());

    ParIds pids2;
    pids2.insert_all(pids);
    BOOST_CHECK(pids2.size()==197);
    i=0;
    for(it=pids2.begin(); it!=pids2.end(); ++it){
        if(i!=10 && i!=102 && i!=198){
            BOOST_CHECK(*it == vpids[i]);
            ++i;
        }
    }
    BOOST_CHECK(pids == pids2);



}

BOOST_AUTO_TEST_SUITE_END()

