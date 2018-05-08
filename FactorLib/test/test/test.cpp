// test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>
#include "mpir.h"
#include "factorlib.h"

int main()
{
	std::vector<unsigned int> primes;
	unsigned long long i;
	mpz_t div1, div2, n, c;
	mpz_init( div1 );
	mpz_init( div2 );
	mpz_init(n);
	mpz_init(c);
	mpz_set_str( n , "4999486012441", 10);
	mpz_set_ui( c, 41 );
	//FactorLib::FactorLib::PollardRho( div1, n, c, 10000 );
	//FactorLib::FactorLib::TonelliShanks( div, n, c );
	FactorLib::FactorLib::QuadraticSieve( div1, div2, n, 127, 10000 );
	i = mpz_get_ui( div1 );
	std::cout << i << std::endl;
	mpz_clear( div1 );
	mpz_clear( div2 );
	mpz_clear( n );
	mpz_clear( c );
	std::cin >> i;
}

