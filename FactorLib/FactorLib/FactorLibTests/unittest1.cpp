#include "stdafx.h"
#include "CppUnitTest.h"
#include <math.h>
#include "../FactorLib/factorlib.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace FactorLibTests
{		
	TEST_CLASS( UnitTest1 )
	{
	public:
		
		TEST_METHOD( TestGCD )
		{
			mpz_t out, a, b;
			mpz_init( out );
			mpz_init( a   ); 
			mpz_init( b   );

			mpz_set_ui( a, 835 );
			mpz_set_ui( b, 4898725 );
			FactorLib::FactorLib::GCD( out, a, b );
			Assert::IsTrue( mpz_cmp_ui( out, 5 ) == 0 );

			mpz_set_ui( a, 835879 );
			mpz_set_ui( b, 4898725 );
			FactorLib::FactorLib::GCD( out, a, b );
			Assert::IsTrue( mpz_cmp_ui( out, 1 ) == 0 );

			mpz_set_str( a, "48797987987984", 10 );
			mpz_set_str( b, "464879878942", 10 );
			FactorLib::FactorLib::GCD( out, a, b );
			Assert::IsTrue( mpz_cmp_ui( out, 2 ) == 0 );

			mpz_set_str( a, "9801000039", 10 );
			mpz_set_str( b, "6534000026", 10 );
			FactorLib::FactorLib::GCD( out, a, b );
			Assert::IsTrue( mpz_cmp_ui( out, 3267000013 ) == 0 );

			mpz_set_str( a, "326700001349898797798949849498", 10 );
			mpz_set_str( b, "65340000264646846816165684867", 10 );
			FactorLib::FactorLib::GCD( out, a, b );
			Assert::IsTrue( mpz_cmp_ui( out, 1 ) == 0 );

			mpz_clear( out );
			mpz_clear( a   );
			mpz_clear( b   ); 
		}

		TEST_METHOD( TestSieveOfE )
		{
			std::vector<unsigned long long> primes = FactorLib::FactorLib::SieveOfE(100);
			Assert::IsTrue( (int)primes.size() == 25 );

		    primes = FactorLib::FactorLib::SieveOfE(1000);
			Assert::IsTrue( (int)primes.size() == 168 );

			primes = FactorLib::FactorLib::SieveOfE(10000);
			Assert::IsTrue( (int)primes.size() == 1229 );

			primes = FactorLib::FactorLib::SieveOfE(100000);
			Assert::IsTrue( (int)primes.size() == 9592	 );

			primes = FactorLib::FactorLib::SieveOfE(1000000);
			Assert::IsTrue( (int)primes.size() == 78498 );
		}

		TEST_METHOD( TestModExpo )
		{
			mpz_t ret, a, b, mod;
			mpz_init( ret );
			mpz_init( a    ); 
			mpz_init( b   );
			mpz_init( mod );

			mpz_set_str( a  , "2"    , 10 );
			mpz_set_str( b  , "10000", 10 );
			mpz_set_str( mod, "13"   , 10 );
			FactorLib::FactorLib::ModExpo( ret, a, b, mod );
			Assert::IsTrue( mpz_cmp_ui( ret, 3 ) == 0 );

			mpz_set_str( a  , "6487"  , 10 );
			mpz_set_str( b  , "489797", 10 );
			mpz_set_str( mod, "56"    , 10 );
			FactorLib::FactorLib::ModExpo( ret, a, b, mod );
			Assert::IsTrue( mpz_cmp_ui( ret, 31 ) == 0 );

			mpz_set_str( a  , "5618494984984"  , 10 );
			mpz_set_str( b  , "156489898987987", 10 );
			mpz_set_str( mod, "2658498489"     , 10 );
			FactorLib::FactorLib::ModExpo( ret, a, b, mod );
			Assert::IsTrue( mpz_cmp_ui( ret, 183734746 ) == 0 );

			mpz_set_str( a  , "561849498498416544989"   , 10 );
			mpz_set_str( b  , "154648989956489898987987", 10 );
			mpz_set_str( mod, "2658498489164984"        , 10 );
			FactorLib::FactorLib::ModExpo( ret, a, b, mod );
			Assert::IsTrue( mpz_cmp_ui( ret, 2093119650212917 ) == 0 );

			mpz_set_str( a  , "2"          , 10 );
			mpz_set_str( b  , "18"         , 10 );
			mpz_set_str( mod, "50000000000", 10 );
			FactorLib::FactorLib::ModExpo( ret, a, b, mod );
			Assert::IsTrue( mpz_cmp_ui( ret, 262144 ) == 0 );

			mpz_clear( ret );
			mpz_clear( a   ); 
			mpz_clear( b   );
			mpz_clear( mod );
		}

		TEST_METHOD( TestPrimeTest )
		{
			bool isPrime;
			isPrime = FactorLib::FactorLib::DeterministicRabinMillerPrimeTest(1);
			Assert::IsTrue( !isPrime );
			isPrime = FactorLib::FactorLib::DeterministicRabinMillerPrimeTest(9576890767);
			Assert::IsTrue( isPrime );
			isPrime = FactorLib::FactorLib::DeterministicRabinMillerPrimeTest(949498448897);
			Assert::IsTrue( !isPrime );
			isPrime = FactorLib::FactorLib::DeterministicRabinMillerPrimeTest(1338);
			Assert::IsTrue( !isPrime );
			isPrime = FactorLib::FactorLib::DeterministicRabinMillerPrimeTest(101);
			Assert::IsTrue( isPrime );
		}

		TEST_METHOD( TestTrialDiv )
		{
			std::map<std::string, unsigned long long, FactorLib::FactorMapCmp> factors;
			FactorLib::FactorLib::TrialDiv( 44879, factors );
			unsigned long long result = 1;
			for( auto it = factors.begin(); it != factors.end(); ++it )
			{
				result *= pow( std::stoull(it->first), it->second );
			}
			Assert::IsTrue( result == 44879 );

			factors.clear();
			FactorLib::FactorLib::TrialDiv( 17, factors );
			result = 1;
			for( auto it = factors.begin(); it != factors.end(); ++it )
			{
				result *= pow( std::stoull(it->first), it->second );
			}
			Assert::IsTrue( result == 17 );

			factors.clear();
			FactorLib::FactorLib::TrialDiv( 1000, factors );
			result = 1;
			for( auto it = factors.begin(); it != factors.end(); ++it )
			{
				result *= pow( std::stoull(it->first), it->second );
			}
			Assert::IsTrue( result == 1000 );

			factors.clear();
			FactorLib::FactorLib::TrialDiv( 789, factors );
			result = 1;
			for( auto it = factors.begin(); it != factors.end(); ++it )
			{
				result *= pow( std::stoull(it->first), it->second );
			}
			Assert::IsTrue( result == 789 );

			factors.clear();
			FactorLib::FactorLib::TrialDiv( 7433, factors );
			result = 1;
			for( auto it = factors.begin(); it != factors.end(); ++it )
			{
				result *= pow( std::stoull(it->first), it->second );
			}
			Assert::IsTrue( result == 7433 );
		}

		TEST_METHOD( TestFermat )
		{
			mpz_t n, a, b;
			mpz_init( n );
			mpz_init( a ); 
			mpz_init( b );

			mpz_set_ui( n, 271420541 );
			FactorLib::FactorLib::Fermat( a, b, n );
			Assert::IsFalse( mpz_cmp( a, n ) == 0 );
			Assert::IsFalse( mpz_cmp( b, n ) == 0 );
			mpz_mul( a, a, b);
			Assert::IsTrue( mpz_cmp( a, n ) == 0 );

			mpz_set_ui( n, 27144361209 );
			FactorLib::FactorLib::Fermat( a, b, n );
			Assert::IsFalse( mpz_cmp( a, n ) == 0 );
			Assert::IsFalse( mpz_cmp( b, n ) == 0 );
			mpz_mul( a, a, b);
			Assert::IsTrue( mpz_cmp( a, n ) == 0 );

			mpz_set_ui( n, 4845050131088271 );
			FactorLib::FactorLib::Fermat( a, b, n );
			Assert::IsFalse( mpz_cmp( a, n ) == 0 );
			Assert::IsFalse( mpz_cmp( b, n ) == 0 );
			mpz_mul( a, a, b);
			Assert::IsTrue( mpz_cmp( a, n ) == 0 );

			mpz_set_ui( n, 76437 );
			FactorLib::FactorLib::Fermat( a, b, n );
			Assert::IsFalse( mpz_cmp( a, n ) == 0 );
			Assert::IsFalse( mpz_cmp( b, n ) == 0 );
			mpz_mul( a, a, b);
			Assert::IsTrue( mpz_cmp( a, n ) == 0 );

			mpz_set_ui( n, 1010426592163 );
			FactorLib::FactorLib::Fermat( a, b, n );
			Assert::IsFalse( mpz_cmp( a, n ) == 0 );
			Assert::IsFalse( mpz_cmp( b, n ) == 0 );
			mpz_mul( a, a, b);
			Assert::IsTrue( mpz_cmp( a, n ) == 0 );

			mpz_clear( n );
			mpz_clear( a ); 
			mpz_clear( b );
		}

		TEST_METHOD( TestPollardRho )
		{
			mpz_t div, n, mod;
			mpz_init( div );
			mpz_init( n   );
			mpz_init( mod );

			mpz_set_ui( n, 168494849847);
			FactorLib::FactorLib::PollardRho( div, n, 10000 );
			Assert::IsFalse( mpz_cmp( div, n ) == 0 );
			mpz_mod( mod, n, div );
			Assert::IsTrue( mpz_cmp_ui( mod, 0 ) == 0 );

			mpz_set_ui( n, 16849);
			FactorLib::FactorLib::PollardRho( div, n, 10000 );
			Assert::IsFalse( mpz_cmp( div, n ) == 0 );
			mpz_mod( mod, n, div );
			Assert::IsTrue( mpz_cmp_ui( mod, 0 ) == 0 );

			mpz_set_ui( n, 168494849847489799);
			FactorLib::FactorLib::PollardRho( div, n, 10000 );
			Assert::IsFalse( mpz_cmp( div, n ) == 0 );
			mpz_mod( mod, n, div );
			Assert::IsTrue( mpz_cmp_ui( mod, 0 ) == 0 );

			mpz_set_ui( n, 1684948498471513);
			FactorLib::FactorLib::PollardRho( div, n, 10000 );
			Assert::IsFalse( mpz_cmp( div, n ) == 0 );
			mpz_mod( mod, n, div );
			Assert::IsTrue( mpz_cmp_ui( mod, 0 ) == 0 );

			mpz_set_ui( n, 168449847);
			FactorLib::FactorLib::PollardRho( div, n, 10000 );
			Assert::IsFalse( mpz_cmp( div, n ) == 0 );
			mpz_mod( mod, n, div );
			Assert::IsTrue( mpz_cmp_ui( mod, 0 ) == 0 );

			mpz_clear( div );
			mpz_clear( n   );
			mpz_clear( mod );
		}

		TEST_METHOD( TestPMinus )
		{
			mpz_t ret, n, mod;
			mpz_init( ret );
			mpz_init( n   );
			mpz_init( mod );

			mpz_set_ui( n, 24748738925);
			FactorLib::FactorLib::Pminus( ret, n );
			Assert::IsFalse( mpz_cmp( ret, n ) == 0 );
			mpz_mod( mod, n, ret );
			Assert::IsTrue( mpz_cmp_ui( mod, 0 ) == 0 );

			mpz_set_ui( n, 2474873);
			FactorLib::FactorLib::Pminus( ret, n );
			Assert::IsFalse( mpz_cmp( ret, n ) == 0 );
			mpz_mod( mod, n, ret );
			Assert::IsTrue( mpz_cmp_ui( mod, 0 ) == 0 );

			mpz_set_ui( n, 247487389254987);
			FactorLib::FactorLib::Pminus( ret, n );
			Assert::IsFalse( mpz_cmp( ret, n ) == 0 );
			mpz_mod( mod, n, ret );
			Assert::IsTrue( mpz_cmp_ui( mod, 0 ) == 0 );

			mpz_set_ui( n, 6744612902657236153);
			FactorLib::FactorLib::Pminus( ret, n );
			Assert::IsTrue( mpz_cmp( ret, n ) == 0 );
			mpz_mod( mod, n, ret );
			Assert::IsTrue( mpz_cmp_ui( mod, 0 ) == 0 );

			mpz_set_str( n, "486933226074853156190591", 10 );
			FactorLib::FactorLib::Pminus( ret, n );
			Assert::IsFalse( mpz_cmp( ret, n ) == 0 );
			mpz_mod( mod, n, ret );
			Assert::IsTrue( mpz_cmp_ui( mod, 0 ) == 0 );

			mpz_clear( ret );
			mpz_clear( n   );
			mpz_clear( mod );
		}

		TEST_METHOD( TestTonelliShanks )
		{
			mpz_t ret, n, mod;
			mpz_init( ret );
			mpz_init( n   );
			mpz_init( mod );

			mpz_set_str( n, "165486878977", 10 );
			mpz_set_ui( mod, 17);
			FactorLib::FactorLib::TonelliShanks( ret, n, mod );
			Assert::IsTrue( mpz_cmp_ui( ret, 8 ) == 0 || mpz_cmp_ui( ret, 9 ) == 0 );

			mpz_set_str( n, "16548687897487977", 10 );
			mpz_set_ui( mod, 71);
			FactorLib::FactorLib::TonelliShanks( ret, n, mod );
			Assert::IsTrue( mpz_cmp_ui( ret, 12 ) == 0 || mpz_cmp_ui( ret, 59 ) == 0 );

			mpz_set_str( n, "14657", 10 );
			mpz_set_ui( mod, 11);
			FactorLib::FactorLib::TonelliShanks( ret, n, mod );
			Assert::IsTrue( mpz_cmp_ui( ret, 4 ) == 0 || mpz_cmp_ui( ret, 7 ) == 0 );

			mpz_set_str( n, "4894849897984894674897", 10 );
			mpz_set_ui( mod, 1500450271);
			FactorLib::FactorLib::TonelliShanks( ret, n, mod );
			Assert::IsTrue( mpz_cmp_ui( ret, 209879152 ) == 0 || mpz_cmp_ui( ret, 1290571119 ) == 0 );

			mpz_set_str( n, "9879789499498484984848948489494494949", 10 );
			mpz_set_ui( mod, 10000000002065383);
			FactorLib::FactorLib::TonelliShanks( ret, n, mod );
			Assert::IsTrue( mpz_cmp_ui( ret, 2206178374379924 ) == 0 || mpz_cmp_ui( ret, 7793821627685459 ) == 0 );

			mpz_clear( ret );
			mpz_clear( n   );
			mpz_clear( mod );
		}

		TEST_METHOD( TestQuadraticSieve )
		{
			mpz_t div1, div2, n;
			mpz_init( div1 );
			mpz_init( div2 );
			mpz_init( n );
			
			mpz_set_str( n, "1164653", 10 );
			FactorLib::FactorLib::QuadraticSieve( div1, div2, n );
			Assert::IsFalse( mpz_cmp( div1, n ) == 0 );
			Assert::IsFalse( mpz_cmp( div2, n ) == 0 );
			mpz_mul( div1, div1, div2 );
			Assert::IsTrue( mpz_cmp( div1, n) == 0 );

			mpz_set_str( n, "151616488497", 10 );
			FactorLib::FactorLib::QuadraticSieve( div1, div2, n );
			Assert::IsFalse( mpz_cmp( div1, n ) == 0 );
			Assert::IsFalse( mpz_cmp( div2, n ) == 0 );
			mpz_mul( div1, div1, div2 );
			Assert::IsTrue( mpz_cmp( div1, n) == 0 );

			mpz_set_str( n, "498484841891891897", 10 );
			FactorLib::FactorLib::QuadraticSieve( div1, div2, n );
			Assert::IsFalse( mpz_cmp( div1, n ) == 0 );
			Assert::IsFalse( mpz_cmp( div2, n ) == 0 );
			mpz_mul( div1, div1, div2 );
			Assert::IsTrue( mpz_cmp( div1, n) == 0 );

			mpz_set_str( n, "484987897984984949898498945", 10 );
			FactorLib::FactorLib::QuadraticSieve( div1, div2, n );
			Assert::IsFalse( mpz_cmp( div1, n ) == 0 );
			Assert::IsFalse( mpz_cmp( div2, n ) == 0 );
			mpz_mul( div1, div1, div2 );
			Assert::IsTrue( mpz_cmp( div1, n) == 0 );

			mpz_set_str( n, "993474687978897714878978987989", 10  );
			FactorLib::FactorLib::QuadraticSieve( div1, div2, n );
			Assert::IsFalse( mpz_cmp( div1, n ) == 0 );
			Assert::IsFalse( mpz_cmp( div2, n ) == 0 );
			mpz_mul( div1, div1, div2 );
			Assert::IsTrue( mpz_cmp( div1, n) == 0 );

			mpz_clear( div1 );
			mpz_clear( div2 );
			mpz_clear( n    );
		}

		TEST_METHOD( TestFactorize )
		{
			mpz_t n, result, div;
			mpz_init( n );
			mpz_init( result );
			mpz_init( div );
			std::map<std::string, unsigned long long, FactorLib::FactorMapCmp> factors;

			mpz_set_ui( result, 1 );
			mpz_set_str( n, "139", 10 );
			FactorLib::FactorLib::Factorize( n, factors );
			for( auto it = factors.begin(); it != factors.end(); ++it )
			{
				mpz_set_str( div, it->first.c_str(), 10 );
				mpz_pow_ui( div, div, it->second );
				mpz_mul( result, result, div );
			}
			Assert::IsTrue( mpz_cmp_ui( result, 139 ) == 0 );

			factors.clear();
			mpz_set_ui( result, 1 );
			mpz_set_str( n, "1164653", 10 );
			FactorLib::FactorLib::Factorize( n, factors );
			for( auto it = factors.begin(); it != factors.end(); ++it )
			{
				mpz_set_str( div, it->first.c_str(), 10 );
				mpz_pow_ui( div, div, it->second );
				mpz_mul( result, result, div );
			}
			Assert::IsTrue( mpz_cmp_ui( result, 1164653 ) == 0 );

			factors.clear();
			mpz_set_ui( result, 1 );
			mpz_set_str( n, "151616488497", 10 );
			FactorLib::FactorLib::Factorize( n, factors );
			for( auto it = factors.begin(); it != factors.end(); ++it )
			{
				mpz_set_str( div, it->first.c_str(), 10 );
				mpz_pow_ui( div, div, it->second );
				mpz_mul( result, result, div );
			}
			Assert::IsTrue( mpz_cmp_ui( result, 151616488497 ) == 0 );

			factors.clear();
			mpz_set_ui( result, 1 );
			mpz_set_str( n, "498484841891891897", 10 );
			FactorLib::FactorLib::Factorize( n, factors );
			for( auto it = factors.begin(); it != factors.end(); ++it )
			{
				mpz_set_str( div, it->first.c_str(), 10 );
				mpz_pow_ui( div, div, it->second );
				mpz_mul( result, result, div );
			}
			Assert::IsTrue( mpz_cmp_ui( result, 498484841891891897 ) == 0 );

			factors.clear();
			mpz_set_ui( result, 1 );
			mpz_set_str( n, "484987897984984949898498945", 10 );
			FactorLib::FactorLib::Factorize( n, factors );
			for( auto it = factors.begin(); it != factors.end(); ++it )
			{
				mpz_set_str( div, it->first.c_str(), 10 );
				mpz_pow_ui( div, div, it->second );
				mpz_mul( result, result, div );
			}
			mpz_set_str( n, "484987897984984949898498945", 10 );
			Assert::IsTrue( mpz_cmp( result, n ) == 0 );

			mpz_clear( result );
			mpz_clear( div );
			mpz_clear( n );
		}

	};
}