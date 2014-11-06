#include <boost/python.hpp>
#include "DmpAlgBgoPID.h"

BOOST_PYTHON_MODULE(libDmpAlgBgoPID){
  using namespace boost::python;

  class_<DmpAlgBgoPID,boost::noncopyable,bases<DmpVAlg> >("DmpAlgBgoPID",init<>());
}
