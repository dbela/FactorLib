#pragma once

#ifndef FACTORLIB_H
#define FACTORLIB_H

#include "..\..\mpir\lib\x64\Release\mpir.h"
#include <vector>
#include <map>
#include <string>
#include <mutex>
#include <thread>

namespace FactorLib
{
	static bool CompFactorMap( std::string a, std::string b )
	{
		if( a.length() != b.length() )
		{
			return a.length() < b.length();
		}
		else
		{
			for( int i = 0; i < a.length(); ++i )
			{
				if( a[i] != b[i] )
				{
					return a[i] < b[i];
				}
			}
		}
		return true;
	}

	struct FactorMapCmp 
	{
		bool operator()(std::string a, std::string b) const
		{
			return CompFactorMap(a, b);
		}
	};

	typedef unsigned long long uLongLong;

	typedef std::map<std::string, uLongLong, FactorMapCmp> factorMap;


	class FactorLib
	{
		public:

			static void                           GCD( mpz_t out, mpz_t a, mpz_t b);
			 							          
			static std::vector<uLongLong>         SieveOfE( uLongLong n );
										          
			static void                           ModExpo( mpz_t ret, mpz_t a, mpz_t b, mpz_t mod );
			
			static uLongLong                      SPRP( uLongLong n, uLongLong b );

			static bool                           DeterministicRabinMillerPrimeTest( uLongLong n );
			
			static void                           TrialDiv( uLongLong n, factorMap &factors );
										          
			static void                           Fermat( mpz_t a, mpz_t b, mpz_t n );
										          
			static void                           PollardRho( mpz_t divisor, mpz_t n, unsigned int maxpolynoms );
										          
			static void                           Pminus( mpz_t ret, mpz_t n );
										          							          
			static bool                           EulerCriterion( mpz_t n, mpz_t mod);
										          
			static uLongLong                      LesserPrimesCount( uLongLong n );
						                          
			static void                           TonelliShanks( mpz_t ret, mpz_t n, mpz_t mod);
										      						          
			static bool                           CanBeFactoredOnBase( std::vector<uLongLong> &vecFactor, std::vector<long> &FactorBase, mpz_t n );
										          
			static double                         GetT( mpz_t n );                 

			static void                           SieveOfQCheckNumber( mpz_t* smoothBases, mpz_t *arrSieve, std::vector<std::vector<uLongLong> > &vecFactors,std::vector<long> &FactorBase, std::vector<float> &vecCheck, mpz_t LowerBound, mpz_t n, double CloseNUF, int start, int add);
			
			static mpz_t*                         SieveOfQ( mpz_t* smoothBases, std::vector<std::vector<uLongLong> > &vecFactors,std::vector<long> &FactorBase, mpz_t n, uLongLong B, uLongLong size );
										          
			static uLongLong                      CheckMatrix( std::vector<std::vector<int> > &ToCheck );
										          
			static void                           AddRowBinary( std::vector<std::vector<int> > &ToElim,  uLongLong row, uLongLong col );
										          
			static uLongLong                      BinaryGaussElimination( std::vector<std::vector<int> > &ToElim );

			static std::vector<std::vector<int> > GetBinaryMatrix( std::vector<std::vector<uLongLong> > vecFactors );

			static void                           GetMultiplier( mpz_t nMultiplier, mpz_t n );
										          
			static void                           GetFactorBase( std::vector<long> &FactorBase, mpz_t n, uLongLong B );

			static void                           FailedSieve( uLongLong &B, uLongLong &size, mpz_t n );

			static uLongLong                      GetB( mpz_t n );
			
			static uLongLong                      GetSize( mpz_t n );
										          
			static void                           QuadraticSieve( mpz_t div1, mpz_t div2, mpz_t n, uLongLong &baseSize, int &GaussNum ,uLongLong B = 0, uLongLong size = 0 );

			static void                           Factorize( mpz_t n, factorMap &factors );
			
			static std::string                    RunFactorize( std::string strNum );

			static std::string                    RunFermat( std::string strNum );

			static std::string                    RunPMinus( std::string strNum );

			static std::string                    RunPollardRho( std::string strNum );

			static std::string                    RunQuadraticSieve( std::string strNum, std::string B = "", std::string size ="" );

			static std::string                    RunPrimeTest( std::string strNum );

			static std::vector<uLongLong>         RunSieveOfE( std::string strNum );
	};

	static std::vector<unsigned long long> vecMillionPrimes = FactorLib::SieveOfE(1000000);

	static std::mutex m;

}


#endif //FACTORLIB_H