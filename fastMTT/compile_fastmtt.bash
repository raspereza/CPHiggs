#!/bin/bash
#c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) fastmtt_cpp.cpp -o ${CMSSW_BASE}/lib/el9_amd64_gcc12/libfastmtt_cpp.so
c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) fastmtt_cpp.cpp -o fastmtt_cpp$(python3-config --extension-suffix)
