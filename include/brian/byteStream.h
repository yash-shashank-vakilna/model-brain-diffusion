#ifndef BYTESTREAM_H
#define BYTESTREAM_H

/*
 *
 * byteStream.h: serialization of data types into and from memory buffers
 * BRIAN Software Package Version 3.0
 *
 * $Id: byteStream.h 416 2016-10-11 01:59:41Z frithjof $
 *
 * 0.10 (17/09/16): initial version
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions and classes for the serialization of data types.
*/

namespace byteStream {
	enum class Endian { Little, Big };
	using BigEndian = std::integral_constant<Endian, Endian::Big>;
	using LittleEndian = std::integral_constant<Endian, Endian::Little>;

	template<typename T, size_t N>
	struct swapEndian {
		void	operator()(T&) { } };
	template<typename T>
	struct swapEndian<T,8> {
		void	operator()(T& x)
			{ union EightBytes { T x; uint8_t y[8]; };
			  EightBytes b; b.x = x;
			  std::swap(b.y[0],b.y[7]); std::swap(b.y[1],b.y[6]);
			  std::swap(b.y[2],b.y[5]); std::swap(b.y[3],b.y[4]);
			  x = b.x; } };
	template<typename T>
	struct swapEndian<T,4> {
		void	operator()(T& x)
			{ union FourBytes { T x; uint8_t y[4]; };
			  FourBytes b; b.x = x;
			  std::swap(b.y[0],b.y[3]); std::swap(b.y[1],b.y[2]);
			  x = b.x; } };
	template<typename T>
	struct swapEndian<T,2> {
		void	operator()(T& x)
			{ union TwoBytes { T x; uint8_t y[2]; };
			  TwoBytes b; b.x = x; std::swap(b.y[0],b.y[1]); x = b.x; } };
	template<typename T>
	void swap(T&, std::true_type)
		{ }
	template<typename T>
	void swap(T& val, std::false_type)
		{ swapEndian<T, sizeof(T)>()(val); }

	template<typename endianType> class istr {
		std::vector<unsigned char> vec;
		size_t	ix;
		endianType et;
	public:
		istr(const unsigned char* m = nullptr, size_t s = 0)
			: vec(), ix{0}, et()
			{ vec.assign(m, m+s); }
		bool	eof() const
			{ return ix >= vec.size(); }
		template<typename T> void read(T& t)
			{ assert(ix+sizeof(T) <= vec.size());
			  memcpy(&t,&vec[ix],sizeof(T));
			  byteStream::swap(t,et); ix += sizeof(T); }
		void	read(unsigned char* p, size_t n)
			{ assert(ix+n <= vec.size());
			  memcpy(p,&vec[ix],n);
			  ix += n; }
		template<typename T> void readArray(T* t, const size_t n, const repnType r)
			{ switch (r) {
			  case repnType::bit:
				for (size_t i = 0; i < n; ) {
					unsigned char c = 0; read(c);
					for (int p = 7; p >= 0 && i < n; p--) 
						t[i++] = (c >> p) & 1; };
				break;
			  case repnType::byte:
				for (size_t i = 0; i < n; i++)
			  		{ unsigned char c = 0; read(c); t[i] = T(c); }
				break;
			  case repnType::shrt:
				for (size_t i = 0; i < n; i++)
			  		{ short int c = 0; read(c); t[i] = T(c); }
				break;
			  case repnType::sint:
				for (size_t i = 0; i < n; i++)
			  		{ int c = 0; read(c); t[i] = T(c); }
				break;
			  case repnType::flt:
				for (size_t i = 0; i < n; i++)
			  		{ float c = 0; read(c); t[i] = T(c); }
				break;
			  default: throw optException("Unsupported representation %d", static_cast<int>(r)); } }
	};

	template<typename endianType> class ostr {
		std::vector<unsigned char> vec;
		endianType et;
	public:
		ostr()
			: vec(), et() {}
		template<typename T> void save(const T& t)
			{ std::vector<unsigned char> v(sizeof(T)); T u{t}; byteStream::swap(u,et);
			  memcpy(&v[0],&u,sizeof(T));
			  vec.insert(vec.end(),v.begin(),v.end()); }
		template<typename T> void copy(const T& t, const size_t ix)
			{ std::vector<unsigned char> v(sizeof(T)); T u{t}; byteStream::swap(u,et);
			  memcpy(&v[0],&u,sizeof(T));
			  std::copy(v.begin(),v.end(),vec.begin()+ix); }
		void	save(const unsigned char* p, size_t n)
			{ for (size_t i = 0; i < n; i++) vec.push_back(p[i]); }
		const std::vector<unsigned char>& data()
			{ return vec; }
		template<typename T> void saveArray(const T* t, const size_t n, repnType& r)
			{ size_t ix = 0;
			  if (r == repnType::unknown) {
				if (typeid(T) != typeid(bool)) {
					vec.resize(n*sizeof(T));
			  		for (unsigned int i = 0; i < n; i++) {
						copy(t[i], ix); ix += sizeof(T); };
					if (typeid(T) == typeid(unsigned char))
						r = repnType::byte;
					else if (typeid(T) == typeid(short int))
						r = repnType::shrt;
					else if (typeid(T) == typeid(unsigned int))
						r = repnType::sint;
					else if (typeid(T) == typeid(float))
						r = repnType::flt;
					return; }
				else r = repnType::bit; };				// bool must be packed below
			  switch (r) {
			  case repnType::bit: {
				unsigned int nb = n/8; if (n%8) nb++; vec.resize(nb);
				for (size_t i = 0; i < n; ) { unsigned char c = 0;
					for (int p = 7; p >= 0 && i < n; i++, p--) 
						if (t[i]) c |= (1 << p);
					vec[ix++] = c; } }
				break;
			  case repnType::byte: {
				vec.resize(n*sizeof(unsigned char));
				for (size_t i = 0; i < n; i++) {
					double v = static_cast<double>(t[i]);
			  		v = std::max(v,static_cast<double>(std::numeric_limits<unsigned char>::min()));
					v = std::min(v,static_cast<double>(std::numeric_limits<unsigned char>::max()));
					vec[ix++] = static_cast<unsigned char>(v); } }
				break;
			  case repnType::shrt: {
				vec.resize(n*sizeof(short int));
				for (size_t i = 0; i < n; i++) {
					double v = static_cast<double>(t[i]);
					v = std::max(v,static_cast<double>(std::numeric_limits<short int>::min()));
					v = std::min(v,static_cast<double>(std::numeric_limits<short int>::max()));
					const short int c = static_cast<short int>(v);
					copy(c, ix); ix += sizeof(short int); } }
				break;
			  case repnType::sint: {
				vec.resize(n*sizeof(int));
				for (size_t i = 0; i < n; i++) {
					double v = static_cast<double>(t[i]);
					v = std::max(v,static_cast<double>(std::numeric_limits<int>::min()));
					v = std::min(v,static_cast<double>(std::numeric_limits<int>::max()));
					const int c = static_cast<int>(v);
					copy(c, ix); ix += sizeof(int); } }
				break;
			  case repnType::flt: {
				vec.resize(n*sizeof(float));
				for (size_t i = 0; i < n; i++) {
					double v = static_cast<double>(t[i]);
					const float c = static_cast<float>(v);
					copy(c, ix); ix += sizeof(float); } }
				break;
			  default: throw optException("Unsupported representation %d", static_cast<int>(r)); } }
	};
};

using endianType = std::is_same<byteStream::BigEndian,byteStream::LittleEndian>;	// set data endianness as first and system endianness as second parameter
using is = byteStream::istr<endianType>;
using os = byteStream::ostr<endianType>;

#endif

