#include "factorlib.h"
#include <math.h>
#include <iostream>

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

		mpz_set_ui( ret, 1 );
		while( mpz_cmp_ui( b, 0 ) != 0 )
		{
			if( mpz_odd_p( b ) != 0 )
			{
				mpz_mul( ret, ret, a );
				mpz_mod( ret, ret, mod );
			}
			mpz_fdiv_q_ui( b, b, 2 );
			mpz_mul( a, a, a );
			mpz_mod( a, a, mod );
		}
	}

	std::vector<uLongLong> FactorLib::SieveOfE( uLongLong n )
	{
		uLongLong size = (uLongLong)ceil( sqrt( n ) );
		bool *array = new bool[n];
		for( uLongLong i = 0; i < n; ++i )
		{
			array[i] = true;
		}
		std::vector<uLongLong> primes;
		
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

	std::vector<uLongLong> FactorLib::TrialDiv( uLongLong n )
	{
		std::vector<uLongLong> factors;
		
		if( n > 1000000 )
			return factors;


		uLongLong d;
		for( d = 2; d <= 3; ++d )
		{
			uLongLong e = 0;
			while( n % d == 0)
			{
				n = n / d;
				e++;
			}
			for( uLongLong i = 0; i < e; ++i )
			{
				factors.push_back(d);
			}
		}

		d = 5;
		uLongLong add = 2;

		while( d*d <= n )
		{
			uLongLong e = 0;
			while( n % d == 0)
			{
				n = n / d;
				e++;
			}
			for( uLongLong i = 0; i < e; ++i )
			{
				factors.push_back(d);
			}

			d += add;
		}
		
		if( n != 1 && d*d > n)
		{
			factors.push_back( n );
		}

		return factors;
	}

	void FactorLib::Fermat( mpz_t a, mpz_t b, mpz_t n )
	{
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

	void FactorLib::PollardRho( mpz_t divisor, mpz_t n, unsigned int max )
	{
		mpz_t x1, x2, product;
		mpz_init( x1 );
		mpz_init( x2 );
		mpz_init( product );
		mpz_set_ui( divisor, 1);
		int c = 1;
		unsigned int primeIndex = 0;
		while( mpz_cmp_ui( divisor, 1) == 0 && primeIndex + 1 < vecMillionPrimes.size() )
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

	void FactorLib::Pminus( mpz_t ret, mpz_t n, mpz_t c, uLongLong max )
	{
		mpz_t pos, tmpGCD;
		mpz_init( pos );
		mpz_init( tmpGCD );
		for( uLongLong i = 0; i < max; ++i )
		{
			mpz_set_ui( pos, i );
			ModExpo( c, c, pos, n );
			if( i % 10 == 0 )
			{
				mpz_sub_ui( tmpGCD, c, 1 );
				GCD( ret, tmpGCD, n );
				if( mpz_cmp_ui( ret, 1 ) > 0 )
				{
					return;
				}
			}
		}

		mpz_clear( pos );
		mpz_clear( tmpGCD );
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

	void FactorLib::SeqV( mpz_t ret, mpz_t h, mpz_t n, uLongLong i )
	{
		if( i = 1 )
		{
			mpz_set( ret, h );
		}
		else if( i = 2 )
		{
			mpz_t tmpMul;
			mpz_init( tmpMul );

			mpz_mul_ui( ret, n, 2 );
			mpz_mul( tmpMul, h, h );
			mpz_sub( ret, tmpMul, ret );
			
			mpz_clear( tmpMul );
		}
		else
		{
			mpz_t tmpRet;
			mpz_init( tmpRet );

			SeqV( tmpRet, h, n, i - 1 );

			mpz_mul( ret, h, tmpRet );

			SeqV( tmpRet, h, n, i - 2 );
			mpz_mul( tmpRet, n, tmpRet );

			mpz_add( ret, ret, tmpRet );

			mpz_clear( tmpRet );
		}
	}

	void FactorLib::CongruenceSolvingWithLegendreSymbol( mpz_t ret, mpz_t n, mpz_t mod )
	{
		mpz_t h, m, x, v, w, tmpMul;
		mpz_init( h );
		mpz_init( x );
		mpz_init( m );
		mpz_init( v );
		mpz_init( w );
		mpz_init( tmpMul );

		uLongLong j = ( mpz_get_ui( mod ) + 1 ) / 2;

		mpz_mul_ui( h, n, 4 );
		mpz_div( h, h, mod );
		mpz_sub_ui( h, h, 1 );
		mpz_sqrt( h, h );

		//SeqV( ret, h, n, i );
		//mpz_mul_ui( ret, ret, i );
		//mpz_mod( ret, ret, mod );

		mpz_set( m, n );
		mpz_set( v, h );
		mpz_mul( w, h, h );
		mpz_mul_ui( tmpMul, m, 2 );
		mpz_sub( w, w, tmpMul );
		mpz_mod( w, w, mod );
		
		std::vector<uLongLong> b;

		while( j > 0 )
		{
			b.push_back( j % 2 );
			j = (uLongLong)floor( j/2 );
		}

		for( uLongLong i = b.size() - 1; i > 0; --i )
		{
			 mpz_mul( x, v, w );
			 mpz_mul( tmpMul, h, m );
			 mpz_sub( x, x, tmpMul );
			 mpz_mod( x, x, mod );

			 mpz_mul( v, v, v );
			 mpz_mul_ui( tmpMul, m, 2 );
			 mpz_sub( v, v, tmpMul );
			 mpz_mod( v, v, mod );

			 mpz_mul( w, v, w );
			 mpz_mul_ui( tmpMul, n, 2 );
			 mpz_mul( tmpMul, tmpMul, m);
			 mpz_sub( w, w, tmpMul );
			 mpz_mod( w, w, mod );

			 mpz_powm_ui( m, m, 2, mod );

			 if( b[i] == 0 )
			 {
				mpz_set( w, x );
			 }
			 else
			 {
				mpz_set( v, x );
				mpz_mul( m, n, m );
				mpz_mod( m, m, mod );
			 }
		}

		mpz_set( ret, v );

		mpz_clear( h );
		mpz_clear( m );
		mpz_clear( x );
		mpz_clear( v );
		mpz_clear( w );
		mpz_clear( tmpMul );
	}

	void FactorLib::TonelliShanks( mpz_t ret, mpz_t n, mpz_t mod) 
	{
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

	void FactorLib::HenselLemma(mpz_t ret, mpz_t prev, mpz_t n, mpz_t mod)
	{
		mpz_t fderivate;
		mpz_init( fderivate );
		mpz_mul_ui( fderivate, prev, 2 );
		mpz_mod( fderivate, fderivate, mod );
		if (mpz_cmp_ui(fderivate, 0) != 0)
		{
			mpz_mul( ret, prev, prev );
			mpz_sub( ret, ret, n );
			mpz_mul( ret, ret, prev );
			mpz_fdiv_q_ui( ret, ret, 2 );
			mpz_sub( ret, prev, ret );
			mpz_mod( ret, ret, mod );
		}

		mpz_clear( fderivate );
	}

	void FactorLib::SolveQuadraticEQ(mpz_t ret, mpz_t prev, mpz_t n, mpz_t mod, int power)
	{
		if (power > 1)
		{
			HenselLemma( ret, prev, n, mod );
		}
		else
		{
			TonelliShanks( ret, n, mod );
		}
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

	mpz_t* FactorLib::SieveOfQ( mpz_t* smoothBases, std::vector<std::vector<uLongLong> > &vecFactors, std::vector<long> &FactorBase, mpz_t n, uLongLong B )
	{
		mpz_t Square, LowerBound;
		mpz_t FxFunction;
		mpz_t QuadraticEq1, QuadraticEq2;
		mpz_init( Square );
		mpz_init( LowerBound );
		mpz_init( FxFunction );
		mpz_init( QuadraticEq1 );
		mpz_init( QuadraticEq2 );
	
		const uLongLong size = 10000000;
		const uLongLong baseSize = FactorBase.size() + 1 + FactorBase.size()/10;
		std::vector<uLongLong> vecFactor;
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
			vecCheck[ Start + 2*i ] += log10(8);
		}

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
				vecCheck[j*primeBase + Start1 - 1] += log10(primeBase);
			}

			for (uLongLong j = 0; j*primeBase + Start2 - 1 < size; ++j)
			{
				vecCheck[j*primeBase + Start2 -1 ] += log10(primeBase);
			}

			mpz_clear(tmpPrime);
			mpz_clear(tmpMod);
		}

		for (uLongLong i = 0; i < size; ++i)
		{
			vecFactor.clear();
			mpz_t x;
			mpz_init( x );
			mpz_add_ui( x, LowerBound, i + 1 );
			mpz_mul( FxFunction, x, x );
			mpz_sub( FxFunction, FxFunction, n );

			if( vecCheck[i] > CloseNUF )
			{
				if ( CanBeFactoredOnBase( vecFactor, FactorBase, FxFunction ) )
				{	
					mpz_set( arrSieve[k], FxFunction );
					mpz_set( smoothBases[k], x );
					vecFactors.push_back(vecFactor);
					k++;
					if (k == baseSize)
					{
						mpz_clear(x);
						break;
					}
				}				
			}
			mpz_clear(x);
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
					ToElim[i][j] += ToElim[row][j];
					ToElim[i][j] =  ToElim[i][j] % 2;
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
		return 0;
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

	void FactorLib::QuadraticSieve( mpz_t div1, mpz_t div2, mpz_t n, uLongLong B )
	{
		mpz_t nMultiplier;
		mpz_init( nMultiplier );

		GetMultiplier( nMultiplier, n );

		std::vector<long> FactorBase;
		
		GetFactorBase( FactorBase, nMultiplier, B );
		
		int baseSize = (int)(FactorBase.size() + 1 + FactorBase.size()/10);

		std::vector<std::vector<uLongLong> > vecFactors;
		std::vector<std::vector<int> > vecFactorsMod2;
		mpz_t* smoothBases   = new mpz_t[ baseSize ];
		for( int i = 0; i < baseSize ; ++i )
		{
			mpz_init( smoothBases[i] );
		}

		mpz_t* smoothNumbers = SieveOfQ( smoothBases, vecFactors, FactorBase, nMultiplier, B );

		vecFactorsMod2 = GetBinaryMatrix( vecFactors );

		mpz_set_ui( div1, 1 );
		while ( ( mpz_cmp_ui(div1, 1) == 0 || mpz_cmp(div1, nMultiplier) == 0 ) && vecFactorsMod2.size() > FactorBase.size() )
		{
			uLongLong index = BinaryGaussElimination( vecFactorsMod2);

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
}
