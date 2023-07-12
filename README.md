Some tools to help with projects involving linear algebra

Main component is the matrix template class  
Type T must have these basic operations defined : +, *, ==, -T (the opposite)  
Operator ~ must be overloded to give conjugate of T. (a+ib --> a-ib for complex numbers)  
A function myTo_string that returns an std::string to represent T must exist for type T to be able to print the matrix.  
A function myZero that return the zero of the type T must be defined  
