#pragma once

#ifndef FACTORLIB_H
#define FACTORLIB_H

#include "mpir.h"
#include <vector>

typedef unsigned long long uLongLong;

namespace FactorLib
{
	class FactorLib
	{
		public:

			static void                   GCD( mpz_t out, mpz_t a, mpz_t b);
			 
			static std::vector<uLongLong> SieveOfE( uLongLong n );

			static std::vector<uLongLong> TrialDiv( uLongLong n );

			static void                   Fermat( mpz_t a, mpz_t b, mpz_t n );

			static void                   PollardRho( mpz_t divisor, mpz_t n, mpz_t c, unsigned int max );

			static void                   Pminus( mpz_t ret, mpz_t n, mpz_t c, uLongLong max );

			static void                   ModExpo( mpz_t ret, mpz_t a, mpz_t b, mpz_t mod );

			static bool                   EulerCriterion( mpz_t n, mpz_t mod);

			static uLongLong              LesserPrimesCount( uLongLong n );

			static void                   SeqV( mpz_t ret, mpz_t h, mpz_t n, uLongLong i );
						                  
			static void                   CongruenceSolvingWithLegendreSymbol( mpz_t ret, mpz_t n, mpz_t mod );
						                  
			static void                   TonelliShanks( mpz_t ret, mpz_t n, mpz_t mod);

			static bool                   CanBeFactoredOnBase( std::vector<uLongLong> &vecFactor, std::vector<uLongLong> &FactorBase, mpz_t n );

			static mpz_t*                 SieveOfQ( mpz_t* smoothBases, std::vector<std::vector<uLongLong> > &vecFactors,std::vector<uLongLong> &FactorBase, mpz_t n, uLongLong B );

			static uLongLong              CheckMatrix( std::vector<std::vector<int> > &ToCheck );

			static void                   AddRowBinary( std::vector<std::vector<int> > &ToElim,  uLongLong row, uLongLong col );
			
			static uLongLong              BinaryGaussElimination( std::vector<std::vector<int> > &ToElim );

			static void                   QuadraticSieve( mpz_t div1, mpz_t div2, mpz_t n, uLongLong B );

	};

	static std::vector<unsigned long long> vecMillionPrimes = FactorLib::SieveOfE(1000000);
}


#endif "FACTORLIB_H"