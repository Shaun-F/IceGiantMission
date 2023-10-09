%module icegiant
%include "typemaps.i"
%{
    #define SWIG_FILE_WITH_INIT
    #include "LISA.hpp"
    #include "Binary.hpp"
    #include "Constants.hpp"
    #include "utils.hpp"
    #include "VariablesParameters.hpp"
    
%}
%include "std_string.i" 

%include "include/LISA.hpp"
%include "include/Binary.hpp"
%include "include/Constants.hpp"
%include "include/utils.hpp"
%include "include/VariablesParameters.hpp"




