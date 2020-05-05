// #ifndef FourDimension_h
// #define FourDimension_h
// 
// class FourDimension {
// public:
//   typedef std::vector<int>::reference reference ;
//   typedef std::vector<int>::const_reference const_reference ;
//   
//   FourDimension() : dims(){}
//   
//   FourDimension(SEXP dims) ;
//   
//   FourDimension( const FourDimension& other ) : dims(other.dims){}
//   FourDimension& operator=( const FourDimension& other ) {
//     if( *this != other )
//       dims = other.dims ;
//     return *this ;
//   }
//   FourDimension(const size_t& n1) : dims(1){
//     dims[0] = static_cast<int>(n1) ;
//   }
//   FourDimension(const size_t& n1, const size_t& n2) : dims(2){
//     dims[0] = static_cast<int>(n1) ;
//     dims[1] = static_cast<int>(n2) ;
//   }
//   FourDimension(const size_t& n1, const size_t& n2, const size_t& n3) : dims(3){
//     dims[0] = static_cast<int>(n1) ;
//     dims[1] = static_cast<int>(n2) ;
//     dims[2] = static_cast<int>(n3) ;
//   }
//   FourDimension(const size_t& n1, const size_t& n2, const size_t& n3, const size_t& n4) : dims(4){
//     dims[0] = static_cast<int>(n1) ;
//     dims[1] = static_cast<int>(n2) ;
//     dims[2] = static_cast<int>(n3) ;
//     dims[3] = static_cast<int>(n4) ;
//   }
//   operator SEXP() const ;
//   
//   inline int size() const {
//     return (int) dims.size() ;
//   }
//   inline R_xlen_t prod() const {
//     return std::accumulate( dims.begin(), dims.end(), static_cast<R_xlen_t>(1), std::multiplies<R_xlen_t>() );
//   }
//   inline reference operator[](int i){
//     if( i < 0 || i>=static_cast<int>(dims.size()) ) throw std::range_error("index out of bounds") ;
//     return dims[i] ;
//   }
//   inline const_reference operator[](int i) const{
//     if( i < 0 || i>=static_cast<int>(dims.size()) ) throw std::range_error("index out of bounds") ;
//     return dims[i] ;
//   }
//   
// private:
//   std::vector<int> dims;
// };
// 
// #endif