#include "factorlib.h"
#include <math.h>
#include <chrono>
#include <iostream>

int iSmoothNumbersArraySize;

namespace FactorLib
{

	void FactorLib::GCD( mpz_t out, mpz_t a, mpz_t b )
	{
		mpz_t tmpA, tmpB, tmpSwap;
		mpz_init( tmpA );
		mpz_init( tmpB );
		mpz_init( tmpSwap );

		mpz_set( tmpB, b );
		mpz_set( tmpA, a );

		while( mpz_cmp_ui( tmpB, 0 ) != 0 )
		{
			mpz_set( tmpSwap, tmpB );
			mpz_mod( tmpB, tmpA, tmpB );
			mpz_set( tmpA, tmpSwap );
		}

		mpz_abs( tmpA, tmpA );
		mpz_set( out, tmpA );

		mpz_clear( tmpA );
		mpz_clear( tmpB );
		mpz_clear( tmpSwap );
	}

	void FactorLib::ModExpo( mpz_t ret, mpz_t a, mpz_t b, mpz_t mod )
	{
		mpz_t base;
		mpz_t power;
		mpz_init( base  );
		mpz_init( power );

		mpz_set( base,  a );
		mpz_set( power, b );
		mpz_set_ui( ret, 1 );
		while( mpz_cmp_ui( power, 0 ) != 0 )
		{
			if( mpz_odd_p( power ) != 0 )
			{
				mpz_mul( ret, ret, base   );
				mpz_mod( ret, ret, mod );
			}
			mpz_fdiv_q_ui( power, power, 2 );
			mpz_mul( base, base, base );
			mpz_mod( base, base, mod  );
		}

		mpz_clear( base  );
		mpz_clear( power );
	}

	std::vector<uLongLong> FactorLib::SieveOfE( uLongLong n )
	{
		std::vector<uLongLong> primes;
		if( n == 1 )
		{
			return primes;
		}
		uLongLong size = (uLongLong)ceil( sqrt( n ) );
		bool *array = new bool[n];
		for( uLongLong i = 0; i < n; ++i )
		{
			array[i] = true;
		}
		
		
		primes.push_back(2);

		for( uLongLong i = 3; i < size; i += 2 )
		{
			if( array[i] == true )
			{
				primes.push_back(i);
				for( uLongLong j = i*i; j < n; j += i )
				{
					array[j] = false;
				}
			}
		}
		int start;
		if( size % 2 == 0)
		{
			start = 1 ;
		}
		else
		{
			start = 0;
		}
		for( uLongLong i = size + start; i < n; i += 2 )
		{
			if( array[i] == true )
			{
				primes.push_back(i);
			}
		}

		return primes;
	}

	uLongLong FactorLib::SPRP( uLongLong n, uLongLong b )
	{
		if( n == 2 || n == 3 ) 
			return n;
		if( ( n & 1 ) == 0 || n == 1 ) 
			return 0;
		uLongLong m = n - 1, nm = n - 1;
        uLongLong k = 0;
        while( ( m & 1 ) == 0 ) 
		{ 
			m >>= 1; 
			k++; 
		}
		mpz_t tx, tb, tm, tn;
		mpz_init( tx );
		mpz_init( tb );
		mpz_init( tm );
		mpz_init( tn );

		mpz_set_ui( tb, b);
		mpz_set_ui( tm, m);
		mpz_set_ui( tn, n);
		ModExpo( tx, tb, tm, tn );

		mpz_clear( tb );
		mpz_clear( tm );

        uLongLong x = mpz_get_ui( tx );

        if( x == 1 || x == n - 1 ) 
			return n;
        while( k )
		{
			mpz_mul( tx, tx, tx );
			mpz_mod( tx, tx, tn );
			x = mpz_get_ui( tx );
			if( x == n - 1 ) 
			{
				mpz_clear( tx );
				mpz_clear( tn );
				return n;
			}
			if( x == 1 ) 
			{
				mpz_clear( tx );
				mpz_clear( tn );
				return 0;
			}
			k--;
        }
		mpz_clear( tx );
		mpz_clear( tn );
        return 0;
	}

	bool FactorLib::DeterministicRabinMillerPrimeTest( uLongLong n )
	{
		uLongLong base[]       = { 2,3,5,7,11,13,17,19,23,29,31 };
        uLongLong baseLength   = sizeof( base ) / sizeof( *base );
        uLongLong primeLimits[]=
		{ 
			2047ul, 1373653ul, 25326001ul, 3215031751ul, 2152302898747ul,
			3474749660383ul, 341550071728321ul, 341550071728321ul,
            825123056546413051ul, 3825123056546413051ul, 3825123056546413051ul 
		};

        uLongLong primeLimitsLength = sizeof( primeLimits ) / sizeof( *primeLimits );
        for( uLongLong i = 0; i < baseLength; i++ )
		{
          if( SPRP( n, base[i] ) == 0 ) 
			  return false;
          if( n < primeLimits[i] ) 
			  return true;
        }
        std::cout << "Not for sure." << std::endl;
        return true;
	}
	
	void FactorLib::TrialDiv( uLongLong n, factorMap &factors )
	{	
		if( n > 1000000 )
			return;


		uLongLong d;
		for( uLongLong i = 0; i < vecMillionPrimes.size(); ++i )
		{
			d = vecMillionPrimes[i];
			
			if( n == 1 || d*d > n )
				break;

			uLongLong e = 0;
			while( n % d == 0)
			{
				n = n / d;
				e++;
			}
			if( e > 0 )
			{
				std::string num = std::to_string( d );
				if( factors.find( num ) != factors.end() )
				{
					factors[num] += e;
				}
				else
				{
					factors[num] = e;
				}
			}				
		}
		
		if( n != 1 && d*d > n)
		{
			std::string num = std::to_string( n );
			if( factors.find( num ) != factors.end() )
			{
				factors[num] += 1;
			}
			else
			{
				factors[num] = 1;
			}
		}
	}

	void FactorLib::Fermat( mpz_t a, mpz_t b, mpz_t n )
	{
		if( mpz_even_p( n ) != 0 )
		{
			mpz_set_ui( a, 2    );
			mpz_div_ui( b, n, 2 );
			return;
		}

		mpz_t sqrt, u, v, r;
		mpz_init( sqrt );
		mpz_init( u );
		mpz_init( v );
		mpz_init( r );
		mpz_sqrt( sqrt, n ),
		mpz_add_ui( sqrt, sqrt, 1 );
		mpz_mul_ui( u, sqrt, 2) ;
		mpz_add_ui( u, u, 1 );
		mpz_set_ui( v, 1 );
		mpz_mul( r, sqrt, sqrt );
		mpz_sub( r, r, n );
		
		while( mpz_cmp_ui( r, 0 ) != 0 )
		{
			if( mpz_cmp_ui( r, 0 ) > 0 )
			{
				while( mpz_cmp_ui( r, 0 ) > 0)
				{
					mpz_sub( r, r, v );
					mpz_add_ui( v, v, 2);
				}
			}
			if( mpz_cmp_ui( r, 0 ) < 0 )
			{
				mpz_add( r, r, u );
				mpz_add_ui( u, u, 2 );
			}
		}

		mpz_add( a, u, v);
		mpz_sub_ui( a, a, 2);
		mpz_div_ui( a, a, 2);

		mpz_sub( b, u, v);
		mpz_div_ui( b, b, 2);
		
		mpz_clear( sqrt );
		mpz_clear( u );
		mpz_clear( v );
		mpz_clear( r );
	}

	void FactorLib::PollardRho( mpz_t divisor, mpz_t n, unsigned int maxpolynoms )
	{
		if( mpz_even_p( n ) != 0 )
		{
			mpz_set_ui( divisor, 2 );
			return;
		}
		unsigned int max = 10000;
		mpz_t x1, x2, product;
		mpz_init( x1 );
		mpz_init( x2 );
		mpz_init( product );
		mpz_set_ui( divisor, 1);
		int c = 1;
		unsigned int primeIndex = 0;
		while( mpz_cmp_ui( divisor, 1) == 0 && primeIndex + 1 < maxpolynoms )
		{
			unsigned int range = 1;
			unsigned int terms = 0;
			mpz_set_ui( x1, 2 );
			mpz_set_ui( x2, 4);
			mpz_add_ui( x2, x2, c );
			mpz_set_ui( product, 1 );
			while( terms <= max )
			{
				for( unsigned int j = 0; j < range; ++j )
				{
					 mpz_powm_ui( x2, x2, 2, n );
					 mpz_add_ui( x2, x2, c );
					 mpz_mod( x2, x2, n);
					 mpz_submul( product, x1, x2 );
					 mpz_mod( product, product, n );
					 terms++;
					 if( terms % 10 == 0 )
					 {
						GCD( divisor, n, product );
						if( mpz_cmp_ui( divisor, 1 ) > 0  )
						{
							mpz_clear( x1 );
							mpz_clear( x2 );
							mpz_clear( product );
							return;
						}
						mpz_set_ui( product, 1 );
					 }
				}
				mpz_set( x1, x2 );
				range = range * 2;
				for( unsigned int j = 0; j < range; ++j )
				{
					mpz_powm_ui( x2, x2, 2, n );
					mpz_add_ui( x2, x2, c );
					mpz_mod( x2, x2, n);
				}
			}
			c = (int)vecMillionPrimes[primeIndex];
			primeIndex++;
		}
		mpz_clear( x1 );
		mpz_clear( x2 );
		mpz_clear( product );
		mpz_set( divisor, n );

	}

	void FactorLib::Pminus( mpz_t ret, mpz_t n )
	{
		if( mpz_even_p( n ) != 0 )
		{
			mpz_set_ui( ret, 2 );
			return;
		}
		unsigned int max = 10000;
		mpz_t pos, c, tmpGCD, exp;
		mpz_init( pos );
		mpz_init( tmpGCD );
		mpz_init( exp );
		mpz_init( c );

		mpz_set_ui( c, 2 );
		mpz_set( exp, c );
		for( uLongLong i = 0; i < max; ++i )
		{
			mpz_set_ui( pos, i+1 );
			ModExpo( exp, c, pos, n );
			if( (i+1) % 10 == 0 )
			{
				mpz_sub_ui( tmpGCD, exp, 1 );
				GCD( ret, tmpGCD, n );
				if( mpz_cmp_ui( ret, 1 ) > 0 )
				{
					return;
				}
			}
			mpz_set( c, exp );
		}

		mpz_clear( exp );
		mpz_clear( pos );
		mpz_clear( tmpGCD );
		mpz_clear( c );
		mpz_set( ret, n );	
	}

	bool FactorLib::EulerCriterion( mpz_t n, mpz_t mod )
	{
		mpz_t tmpNum;
		mpz_init( tmpNum ); 

		mpz_sub_ui( tmpNum, mod, 1 );
		mpz_div_ui( tmpNum, tmpNum, 2 );
		mpz_powm( tmpNum, n, tmpNum, mod );
		
		bool ret;
		if( mpz_cmp_ui( tmpNum, 1 ) == 0 )
		{
			ret = true;
		}
		else
		{
			ret = false;
		}
		mpz_clear( tmpNum );

		return ret;
	}

	uLongLong FactorLib::LesserPrimesCount( uLongLong n )
	{
		if( n < 1000000 )
		{
			int i = 0;
			while( n > vecMillionPrimes[i] )
			{
				++i;
			}
			return i;
		}
		else
		{
			std::vector<uLongLong> primes = SieveOfE( n );
			return primes.size();
		}
	}

	void FactorLib::TonelliShanks( mpz_t ret, mpz_t n, mpz_t mod) 
	{
		if( mpz_cmp_ui( mod , 3) == 0 )
		{
			mpz_set_ui( ret, 1 );
		}

		mpz_t S, Q, z, t;
		mpz_init( S );
		mpz_init( Q );
		mpz_init( z );
		mpz_init( t );

		mpz_sub_ui( Q, mod, 1 );
		while( mpz_even_p( Q ) != 0 )
		{
			mpz_div_ui( Q, Q, 2 );
			mpz_add_ui( S, S, 1 );
		}		
		
		mpz_set_ui( z, 3 );
		while( EulerCriterion( z , mod ) == true && mpz_cmp( z, mod ) < 0 )
		{
			mpz_add_ui( z, z, 1 );
		}

		if( mpz_cmp( z, mod ) >= 0 )
			return;

		mpz_powm( z, z, Q, mod );
		mpz_powm( t, n, Q, mod );
		mpz_add_ui( Q, Q, 1 );
		mpz_div_ui( Q, Q, 2 );
		mpz_powm( ret, n, Q, mod );

		while( true )
		{
			if( mpz_cmp_ui( t, 0 ) == 0 )
			{
				mpz_set_ui( ret, 0 );
				break;
			}

			if( mpz_cmp_ui( t, 1 ) == 0 )
			{
				break;
			}


			uLongLong i = 1;
			mpz_t tmpT;
			mpz_init( tmpT );
			while( mpz_cmp_ui( S, i ) > 0 )
			{
				uLongLong power = (uLongLong)pow(2, i);
				mpz_powm_ui( tmpT, t, power, mod );
				if( mpz_cmp_ui( tmpT, 1 ) == 0 )
					break;
				++i;
			}

			mpz_clear( tmpT );

			if( mpz_cmp_ui( S, i ) == 0 )
			{
				mpz_set_si( ret, -1 );
				break;
			}

			mpz_t b;
			mpz_init( b );
			
			uLongLong power = mpz_get_ui( S );
			power = power - i - 1;
			power = (uLongLong)pow( 2, power);
			mpz_powm_ui( b, z, power, mod );
			mpz_set_ui( S, i );
			mpz_mul( z, b, b );
			mpz_mod( z, z, mod );
			mpz_mul( t, t, z );
			mpz_mod( t, t, mod );
			mpz_mul( ret, ret, b );
			mpz_mod( ret, ret, mod );

			mpz_clear( b );
		}
		mpz_clear( S );
		mpz_clear( Q );
		mpz_clear( z );
		mpz_clear( t );
	}

	bool FactorLib::CanBeFactoredOnBase( std::vector<uLongLong> &vecFactor, std::vector<long> &FactorBase, mpz_t n )
	{
		mpz_t tmp, locN;
		mpz_init( tmp );
		mpz_init( locN );

		mpz_set( locN, n );

		if( mpz_cmp_ui( n, 0 ) < 0 )
		{
			vecFactor.push_back( 1 );
			mpz_mul_si( locN, n, -1 );
		}
		else
		{
			vecFactor.push_back( 0 );
		}

		for( int i = 1; i < FactorBase.size(); ++i )
		{
			uLongLong power = 0;
			uLongLong prime;
			if ( FactorBase[i] == 8 )
			{
				prime = 2;
			}
			else
			{
				prime = FactorBase[i];
			}
			mpz_mod_ui( tmp, locN, prime );
			while( mpz_cmp_ui( tmp, 0 ) == 0 )
			{
				mpz_div_ui( locN, locN, prime );
				mpz_mod_ui( tmp, locN, prime );
				power++;
			}
			vecFactor.push_back( power );
		}

		mpz_clear( tmp );

		if( mpz_cmp_ui( locN, 1 ) == 0 )
		{
			mpz_clear( locN );
			return true;
		}
		else
		{
			mpz_clear( locN );
			return false;
		}
	}

	double FactorLib::GetT(mpz_t n)
	{
		int size = (int)mpz_sizeinbase( n, 10 );
		if (size < 30)
		{
			return 1.5;
		}
		else if (size < 45)
		{
			return 2;
		}
		else
		{
			return 2.6;
		}
	}

	void FactorLib::SieveOfQCheckNumber(mpz_t* smoothBases, mpz_t *arrSieve, std::vector<std::vector<uLongLong> > &vecFactors, std::vector<long> &FactorBase, std::vector<float> &vecCheck, mpz_t LowerBound, mpz_t n, double CloseNUF, int start, int add )
	{
			std::vector<uLongLong> vecFactor;
			mpz_t x;
			mpz_t FxFunction;
			mpz_init( x );
			mpz_init( FxFunction );
			uLongLong baseSize = FactorBase.size() + 1 + FactorBase.size()/10;

			for( int i = start; i < vecCheck.size(); i += add )
			{
				vecFactor.clear();
				if (iSmoothNumbersArraySize >= baseSize)
				{
					mpz_clear(x);
					mpz_clear( FxFunction );
					return;
				}
				if( vecCheck[i] > CloseNUF )
				{
					mpz_add_ui( x, LowerBound, i + 1 );
					mpz_mul( FxFunction, x, x );
					mpz_sub( FxFunction, FxFunction, n );
					if ( CanBeFactoredOnBase( vecFactor, FactorBase, FxFunction ) )
					{	
						m.lock();
						if (iSmoothNumbersArraySize >= baseSize)
						{
							m.unlock();
							mpz_clear(x);
							mpz_clear( FxFunction );
							return;
						}
						mpz_set( arrSieve[iSmoothNumbersArraySize], FxFunction );
						mpz_set( smoothBases[iSmoothNumbersArraySize], x );
						vecFactors.push_back(vecFactor);
						iSmoothNumbersArraySize++;
						m.unlock();
					}
				}
			}
			mpz_clear(x);
			mpz_clear( FxFunction );
		
	}

	mpz_t* FactorLib::SieveOfQ( mpz_t* smoothBases, std::vector<std::vector<uLongLong> > &vecFactors, std::vector<long> &FactorBase, mpz_t n, uLongLong B,  uLongLong size )
	{
		mpz_t Square, LowerBound;
		mpz_t FxFunction;
		mpz_t QuadraticEq1, QuadraticEq2;
		mpz_init( Square );
		mpz_init( LowerBound );
		mpz_init( FxFunction );
		mpz_init( QuadraticEq1 );
		mpz_init( QuadraticEq2 );
	
		const uLongLong baseSize = FactorBase.size() + 1 + FactorBase.size()/10;
		
		std::vector<float> vecCheck(size);
		mpz_t *arrSieve = new mpz_t[baseSize];

		double Target   = (mpz_sizeinbase( n, 10 ) - 1 )/2 + log10(size);
		double CloseNUF = Target - GetT( n )*log10(FactorBase[FactorBase.size()-1]);

		for( int i = 0; i < baseSize; ++i )
		{
			mpz_init( arrSieve[i] );
		}

		for( uLongLong i = 0; i < size; ++i )
		{
			vecCheck[i] = 0;
		}

		int k = 0;	
		mpz_sqrt  ( Square, n );
		mpz_sub_ui( LowerBound, Square, size/2 );

		int Start;
		if( mpz_odd_p(LowerBound) != 0)
		{
			Start = 0;
		}
		else
		{
			Start = 1;
		}
		for( int i = 0; Start + 2 * i < size; ++i )
		{
			vecCheck[ Start + 2*i ] += (float)log10(8);
		}

		long long dur = 0;

		for( uLongLong i = 2; i < FactorBase.size(); ++i )
		{
			mpz_t tmpPrime, tmpMod;
			mpz_init(tmpPrime);
			mpz_init(tmpMod);

			uLongLong primeBase = FactorBase[i] ;
			mpz_set_ui(tmpPrime, primeBase);
			if (primeBase % 2 == 0)
			{
				mpz_set_ui( QuadraticEq1, 1 );
			}
			else
			{
				TonelliShanks( QuadraticEq1, n, tmpPrime );
				mpz_sub( QuadraticEq2, tmpPrime, QuadraticEq1 );
			}
			
			uLongLong QE1 = mpz_get_ui(QuadraticEq1);
			uLongLong QE2 = mpz_get_ui(QuadraticEq2);

			mpz_mod( tmpMod, LowerBound, tmpPrime );
			mpz_sub( tmpMod, QuadraticEq1, tmpMod );
			mpz_mod( tmpMod, tmpMod, tmpPrime );
			uLongLong Start1 = mpz_get_ui( tmpMod );

			mpz_mod( tmpMod, LowerBound, tmpPrime );
			mpz_sub( tmpMod, QuadraticEq2, tmpMod );
			mpz_mod( tmpMod, tmpMod, tmpPrime );
			uLongLong Start2 = mpz_get_ui( tmpMod );
				
			for (uLongLong j = 0; j*primeBase + Start1 - 1< size; ++j)
			{
				vecCheck[j*primeBase + Start1 - 1] += (float)log10(primeBase);
			}

			for (uLongLong j = 0; j*primeBase + Start2 - 1 < size; ++j)
			{
				vecCheck[j*primeBase + Start2 -1 ] += (float)log10(primeBase);
			}

			mpz_clear(tmpPrime);
			mpz_clear(tmpMod);
		}

		iSmoothNumbersArraySize = 0;

		int iAvailableCores = std::thread::hardware_concurrency() - 1;
		std::vector<std::thread> threads(iAvailableCores);

		for( int i = 0; i < iAvailableCores; ++i )
		{
			threads[i] = std::thread( SieveOfQCheckNumber, smoothBases, arrSieve, std::ref(vecFactors), FactorBase, std::ref(vecCheck), LowerBound, n, CloseNUF, i, iAvailableCores );
		}

		for( int i = 0; i < iAvailableCores; ++i )
		{
			threads[i].join();
		}

		
	
		
		mpz_clear( Square );
		mpz_clear( FxFunction );
		
		return arrSieve;
	}

	void FactorLib::AddRowBinary( std::vector<std::vector<int> > &ToElim, uLongLong row, uLongLong col )
	{
		for( uLongLong i = row + 1; i < ToElim.size(); ++i )
		{
			if( ToElim[i][col] == 1 )
			{
				for( uLongLong j = col; j < ToElim[i].size(); ++j )
				{
					ToElim[i][j]^=ToElim[row][j];
				}
			}
		}
	}

	uLongLong FactorLib::CheckMatrix( std::vector<std::vector<int> > &ToCheck )
	{
		for( uLongLong i = 0; i < ToCheck.size(); ++i )
		{
			bool isAllZero = true;
			for( uLongLong j = 0; j < ToCheck.size() - 1 ; ++j )
			{
				isAllZero &= ( ToCheck[i][j] == 0  );
			}
			if( isAllZero == true )
			{
				return i;
			}
		}
		return ToCheck.size();
	}

	uLongLong FactorLib::BinaryGaussElimination( std::vector<std::vector<int> > &ToElim )
	{
		uLongLong result;
		for( uLongLong i = 0; i < ToElim[0].size(); ++i )
		{
			for( uLongLong j = i; j < ToElim.size(); ++j  )
			{
				if( ToElim[j][i] == 1 )
				{
					AddRowBinary( ToElim, j, i );
					if( j != i )
					{
						std::vector<int> tmp = ToElim[j];
						ToElim[j] = ToElim[i];
						ToElim[i] = tmp;
					}
					result = CheckMatrix( ToElim );
					if ( result < ToElim.size() )
					{
						return result;
					}
					break;
				}
			}
		}
		return ToElim.size();
	}

	std::vector<std::vector<int> > FactorLib::GetBinaryMatrix(std::vector<std::vector<uLongLong> > vecFactors)
	{
		std::vector<std::vector<int> > vecFactorsMod2;
		vecFactorsMod2.resize( vecFactors.size() );
		for( uLongLong i = 0; i < vecFactors.size(); ++i )
		{
			for( uLongLong j = 0; j < vecFactors[i].size(); ++j )
			{
				vecFactorsMod2[i].push_back( vecFactors[i][j] % 2 );
			}
			for( uLongLong j = 0; j < vecFactors.size(); ++j )
			{
				if( i == j )
					vecFactorsMod2[i].push_back( 1 );
				else
					vecFactorsMod2[i].push_back( 0 );
			}
		}

		return vecFactorsMod2;
	}

	void FactorLib::GetMultiplier(mpz_t nMultiplier, mpz_t n)
	{
		mpz_t mod8;
		mpz_init( mod8 );

		mpz_mod_ui( mod8, n, 8 );

		if( mpz_cmp_ui(mod8, 3) == 0  )
		{
			mpz_mul_ui( nMultiplier, n, 5 );
		}
		else if (mpz_cmp_ui(mod8, 7) == 0)
		{
			mpz_mul_ui( nMultiplier, n, 7 );
		}
		else if (mpz_cmp_ui(mod8, 5) == 0)
		{
			mpz_mul_ui( nMultiplier, n, 3 );
		}
		else
		{
			mpz_set( nMultiplier, n );
		}

		mpz_clear( mod8 );
	}

	void FactorLib::GetFactorBase(std::vector<long> &FactorBase, mpz_t n, uLongLong B)
	{
		FactorBase.push_back( -1 );
		FactorBase.push_back(  8 );

		uLongLong i = 1;
		uLongLong vecSize = LesserPrimesCount( B );
		while( FactorBase.size() != vecSize && i < vecMillionPrimes.size() )
		{
			mpz_t tmpPrime, tmpBase, tmpBase2;
			mpz_init( tmpPrime );
			mpz_init( tmpBase );
			mpz_init( tmpBase2 );
			mpz_set_ui( tmpPrime, vecMillionPrimes[i] );
			if( EulerCriterion( n, tmpPrime) )
			{
				TonelliShanks( tmpBase, n, tmpPrime );
				mpz_sub( tmpBase2, tmpBase, tmpPrime );

				if( mpz_cmp_ui(tmpBase, B) <= 0 || mpz_cmp_ui( tmpBase2, B ) <= 0 )
					FactorBase.push_back( (long)vecMillionPrimes[i] );
			}
			mpz_clear( tmpPrime  );
			mpz_clear( tmpBase   );
			mpz_clear( tmpBase2  );
			++i;
		}
	}

	uLongLong FactorLib::GetB( mpz_t n )
	{
		size_t size = mpz_sizeinbase( n, 10 );

		if( size <= 10 )
		{
			return 500;
		}
		else if( size <= 15 )
		{
			return 800;
		}
		else if( size <= 20 )
		{
			return 800;
		}
		else if( size <= 25 )
		{
			return 1000;
		}
		else if( size <= 30 )
		{
			return 2000;
		}
		else if( size <= 35 )
		{
			return 5000;
		}
		else if( size <= 40 )
		{
			return 10000;
		}
		else
		{
			return 10000;
		}
		return 0;
	}

	uLongLong FactorLib::GetSize( mpz_t n )
	{
		size_t size = mpz_sizeinbase( n, 10 );
		
		if( size <= 10 )
		{
			return 1000;
		}
		else if( size <= 15 )
		{
			return 10000;
		}
		else if( size <= 20 )
		{
			return 100000;
		}
		else if( size <= 25 )
		{
			return 1000000;
		}
		else if( size <= 30 )
		{
			return 100000000;
		}
		else if( size <= 35 )
		{
			return 300000000;
		}
		else 
		{
			return 700000000;
		}
	}

	void FactorLib::FailedSieve( uLongLong &B, uLongLong& size, mpz_t n )
	{
		size_t numsize = mpz_sizeinbase( n, 10 );

		if( numsize >= 30 )
		{
			B += GetB( n );
		}
		else
		{
			size += 2*size;
		}
	}

	void FactorLib::QuadraticSieve( mpz_t div1, mpz_t div2, mpz_t n, uLongLong &baseSize, int &GaussNum , uLongLong B,  uLongLong size )
	{
		if( mpz_sizeinbase( n, 10 ) < 5 )
			return;

		if( mpz_even_p( n ) != 0 )
		{
			mpz_set_ui( div1, 2    );
			mpz_div_ui( div2, n, 2 );
			return;
		}

		mpz_t nMultiplier;
		mpz_init( nMultiplier );

		GetMultiplier( nMultiplier, n );

		std::vector<long> FactorBase;

		if( B == 0 || size == 0 )
		{
			size = GetSize( n );
			B = GetB( n );
		}

		baseSize = 1;

		std::vector<std::vector<uLongLong> > vecFactors;
		std::vector<std::vector<int> > vecFactorsMod2;
		mpz_t* smoothBases;
		
		mpz_t* smoothNumbers;
		while( vecFactors.size() < baseSize )
		{
			vecFactors.clear();
			FactorBase.clear();
		
			GetFactorBase( FactorBase, nMultiplier, B );

			baseSize = (int)(FactorBase.size() + 1 + FactorBase.size()/10);

			smoothBases = new mpz_t[ baseSize ];

			for( int i = 0; i < baseSize ; ++i )
			{
				mpz_init( smoothBases[i] );
			}

			std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

			smoothNumbers = SieveOfQ( smoothBases, vecFactors, FactorBase, nMultiplier, B , size );

			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
			std::cout << "Sieve" << duration << std::endl;

			FailedSieve( B, size, n );

			//if( B > GetB( n ) * 5 || size > GetSize( n ) * 5 )
			//	break;
		}

		vecFactorsMod2 = GetBinaryMatrix( vecFactors );

		mpz_set_ui( div1, 1 );
		mpz_set( div2, n );
		GaussNum = 0;
		while ( ( mpz_cmp_ui(div1, 1) == 0 || mpz_cmp(div1, n) == 0 ) && vecFactorsMod2.size() > FactorBase.size() )
		{
			GaussNum++;
			std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();

			uLongLong index = BinaryGaussElimination( vecFactorsMod2);

			std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
			auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count();
			std::cout << "Gauss" << duration2 << std::endl;

			mpz_t x, y, tmp;
			mpz_init( y );
			mpz_init( x );
			mpz_init( tmp );

			mpz_set_ui( y, 1 );
			mpz_set_ui( x, 1 );
		

			for( uLongLong i = FactorBase.size(); i < vecFactorsMod2[index].size(); ++i )
			{
				if( vecFactorsMod2[index][i] == 1 )
				{
					mpz_mul( y, y, smoothNumbers[ i - FactorBase.size()] );
					mpz_mul( x, x, smoothBases[ i - FactorBase.size() ] );
				}
			}
			mpz_sqrt( y , y );
			mpz_mod( y, y, nMultiplier );

			mpz_set( tmp, y );
			mpz_sub( y, x, y );
			mpz_add( x, x, tmp);

			GCD( div1, x, n );
			GCD( div2, y, n );

			vecFactorsMod2.erase( vecFactorsMod2.begin() + index );

			mpz_clear( y );
			mpz_clear( x );
			mpz_clear( tmp );
		}
		
		
		for( int i = 0; i < FactorBase.size() + 1 ; ++i  )
		{
			mpz_clear( smoothBases[i] );
			mpz_clear( smoothNumbers[i] );
		}

		
		mpz_clear( nMultiplier );
	}

	void FactorLib::Factorize( mpz_t n, factorMap &factors )
	{
		for( int i = 0; i < vecMillionPrimes.size(); ++i  )
		{
			uLongLong prime = vecMillionPrimes[i];
			int pow2 = 0;
			mpz_t mod;
			mpz_init( mod );
			mpz_mod_ui( mod, n, prime );
			while( mpz_cmp_ui( mod, 0 ) == 0 )
			{	
				mpz_div_ui( n, n, prime );
				pow2++;
				mpz_mod_ui( mod, n, prime );
			}
			if( pow2 > 0 )
			{
				factors[std::to_string(prime)] = pow2;
			}
			if( mpz_cmp_ui( n, 1 ) == 0 )
				return;
		}
		

		int size = (int)mpz_sizeinbase( n, 10 );

		if( size > 19 )
		{	

			mpz_t div1; 
			mpz_t div2; 
			mpz_init( div1 );
			mpz_init( div2 );

			PollardRho( div1, n, 1000 );
			if( mpz_cmp( div1, n ) != 0 )
			{
				mpz_div( div2, n, div1 );

				Factorize( div1, factors );
				Factorize( div2, factors );
				mpz_clear( div1 );
				mpz_clear( div2 );
				return;
			}
			uLongLong basesize;
			int gaussnum;
			QuadraticSieve( div1, div2, n, basesize, gaussnum, 2000, 100000000 );

			if( mpz_cmp( div1, n ) == 0 )
			{
				char* chDiv = new char[mpz_sizeinbase( n, 10 )];
				mpz_get_str( chDiv, 10, n );
				std::string strDiv(chDiv);
				
				if( factors.find( strDiv ) != factors.end() )
				{
					factors[ strDiv ] += 1;
				}
				else
				{
					factors[ strDiv ] = 1;
				}
				mpz_clear( div1 );
				mpz_clear( div2 );
			}
			Factorize( div1, factors );
			Factorize( div2, factors );

			mpz_clear( div1 );
			mpz_clear( div2 );

		}
		else
		{
			uLongLong ullN = mpz_get_ui( n );
			if( DeterministicRabinMillerPrimeTest( ullN ) )
			{
				std::string strDiv = std::to_string( ullN );
				
				if( factors.find( strDiv ) != factors.end() )
				{
					factors[ strDiv ] += 1;
				}
				else
				{
					factors[ strDiv ] = 1;
				}
			}
			else
			{
				if( ullN > 10000 )
				{
					mpz_t div1; 
					mpz_t div2; 
					mpz_init( div1 );
					mpz_init( div2 );
					
					PollardRho( div1, n, 500);
					mpz_div( div2, n, div1 );
					uLongLong baseSize = 0;
					int GaussNum;

					if( mpz_cmp_ui( div1, 1 ) == 0 || mpz_cmp( div1, n ) == 0 )
					{
						QuadraticSieve( div1, div2, n, baseSize, GaussNum );
					}

					Factorize( div1, factors );
					Factorize( div2, factors );
					
					mpz_clear( div1 );
					mpz_clear( div2 );
				}
				else
				{
					TrialDiv( ullN, factors );
				}
			}
		}
	}

	std::string FactorLib::RunFermat( std::string strNum )
	{
		
		mpz_t n, div1, div2;
		mpz_init( n    );
		mpz_init( div1 );
		mpz_init( div2 );

		mpz_set_str( n, strNum.c_str(), 10 );
		
		Fermat( div1, div2, n );

		char* chDiv1 = new char[mpz_sizeinbase( div1, 10 )];
		char* chDiv2 = new char[mpz_sizeinbase( div2, 10 )];
		mpz_get_str( chDiv1, 10, div1 );
		mpz_get_str( chDiv2, 10, div2 );

		mpz_clear( n    );
		mpz_clear( div1 );
		mpz_clear( div2 );
		if( chDiv1 == "1" || chDiv2 == "1" )
		{
			return "A(z) " + strNum + " szám prím.";
		}
		else
		{
			return strNum + " = " + std::string( chDiv1 ) + " * " + std::string( chDiv2 );
		}
		
	}

	std::string FactorLib::RunPMinus( std::string strNum )
	{
		std::string resMsg, strDiv1, strDiv2;
		mpz_t n, div1, div2;
		mpz_init( n    );
		mpz_init( div1 );
		mpz_init( div2 );

		mpz_set_str( n, strNum.c_str(), 10 );

		Pminus( div1, n );

		char* chDiv1 = new char[mpz_sizeinbase( div1, 10 )];
		char* chDiv2 = new char[mpz_sizeinbase( div2, 10 )];
		mpz_get_str( chDiv1, 10, div1 );
		mpz_get_str( chDiv2, 10, div2 );

		mpz_clear( n    );
		mpz_clear( div1 );
		mpz_clear( div2 );
		if( chDiv1 == "1" || chDiv2 == "1" )
		{
			return "A(z) " + strNum + " szám prím.";
		}
		else
		{
			return strNum + " = " + std::string( chDiv1 ) + " * " + std::string( chDiv2 );
		}
	}

	std::string FactorLib::RunPollardRho( std::string strNum )
	{
		std::string resMsg, strDiv1, strDiv2;
		mpz_t n, div1, div2;
		mpz_init( n    );
		mpz_init( div1 );
		mpz_init( div2 );

		mpz_set_str( n, strNum.c_str(), 10 );

		PollardRho( div1, n, 1000 );
		mpz_div( div2, n, div1 );

		char* chDiv1 = new char[mpz_sizeinbase( div1, 10 )];
		char* chDiv2 = new char[mpz_sizeinbase( div2, 10 )];
		mpz_get_str( chDiv1, 10, div1 );
		mpz_get_str( chDiv2, 10, div2 );

		mpz_clear( n    );
		mpz_clear( div1 );
		mpz_clear( div2 );
		if( chDiv1 == "1" || chDiv2 == "1" )
		{
			return "A(z) " + strNum + " szám prím.";
		}
		else
		{
			return strNum + " = " + std::string( chDiv1 ) + " * " + std::string( chDiv2 );
		}
	}
	
	std::string FactorLib::RunQuadraticSieve( std::string strNum, std::string B, std::string size )
	{
		std::string resMsg, strDiv1, strDiv2;
		mpz_t n, div1, div2;
		mpz_init( n    );
		mpz_init( div1 );
		mpz_init( div2 );

		mpz_set_str( n, strNum.c_str(), 10 );

		uLongLong baseSize = 0;
		int GaussNum = 0;
		if( B != "" && size != "" )
		{
			uLongLong base = std::stoull( B );
			QuadraticSieve( div1, div2, n, baseSize, GaussNum, base, std::stoull(size) );
		}
		else
		{
			QuadraticSieve( div1, div2, n, baseSize, GaussNum );
		}
		char* chDiv1 = new char[mpz_sizeinbase( div1, 10 )];
		char* chDiv2 = new char[mpz_sizeinbase( div2, 10 )];
		mpz_get_str( chDiv1, 10, div1 );
		mpz_get_str( chDiv2, 10, div2 );
		
		mpz_clear( n    );
		mpz_clear( div1 );
		mpz_clear( div2 );

		std::string data = "A talált relációk száma: " + std::to_string( baseSize ) + ". A végrehajtott Gauss-eliminációk száma: " + std::to_string( GaussNum ) + ".";

		if( chDiv1 == "1" || chDiv2 == "1" )
		{
			return "A(z) " + strNum + " szám prím.";
		}
		else
		{
			return strNum + " = " + std::string( chDiv1 ) + " * " + std::string( chDiv2 ) + " . " + data;
		}
	}

	std::string FactorLib::RunPrimeTest( std::string strNum )
	{
		mpz_t( n );
		mpz_init( n );
		mpz_set_str( n, strNum.c_str(), 10 );
		bool badNum = (mpz_cmp_ui( n, 18446744073709551615 ) > 0 );
		mpz_clear( n );
		if( badNum )
		{
			return "A szám meghaladja a megadott határértéket.";
		}
		else
		{
			uLongLong n = std::stoull( strNum );
			if( DeterministicRabinMillerPrimeTest( n ) )
			{
				return "A(z) " + strNum + " szám prím.";
			}
			else
			{
				return "A(z) " + strNum + " szám nem prím.";
			}
		}
	}

	std::string FactorLib::RunFactorize( std::string strNum )
	{
		factorMap factors;

		std::string resultMsg = strNum + " = ";
		mpz_t n;
		mpz_init( n );

		mpz_set_str( n, strNum.c_str(), 10 );

		Factorize( n, factors );

		for( auto it = factors.begin(); it != factors.end(); ++it )
		{
			if( it != factors.begin() )
				resultMsg += " * ";
			std::string power = std::to_string( it->second );
			if( power != "1" )
			{
				resultMsg += it->first + "^" + power;
			}
			else
			{
				resultMsg += it->first;
			}
		}
		mpz_clear( n );
		return resultMsg;
	}

	std::vector<uLongLong> FactorLib::RunSieveOfE( std::string strNum )
	{
		mpz_t( n );
		mpz_init( n );
		mpz_set_str( n, strNum.c_str(), 10 );
		bool badNum = (mpz_cmp_ui( n, 18446744073709551615 ) > 0 );
		mpz_clear( n );

		if( badNum )
		{
			std::vector<uLongLong> primes;
			return primes;
		}
		else
		{
			uLongLong n = std::stoull( strNum );
			return SieveOfE( n );
		}
	}
}
