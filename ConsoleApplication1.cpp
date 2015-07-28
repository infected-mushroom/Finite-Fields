// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <set>

struct Isomorphism {
	std::vector<int> preimage;
	std::vector<int> image;
};

int binpow(int a, int n) {
	int res = 1;
	while (n) {
		if (n & 1)
			res *= a;
		a *= a;
		n >>= 1;
	}
	return res;
}


std::vector<std::vector<int> > make_polynomials(int degree, int mod) {
	int q = binpow(mod, degree);
	std::vector<std::vector<int> > result((mod - 1)*q, std::vector<int>(degree + 1, 0));

	for (size_t i = q; i != (q*(mod)); ++i) {
		int num = i;
		int cur_position = degree;
		while (num != 0) {
			int remainder = num % mod;
			result[i - q][cur_position] = remainder;
			num = num / mod;
			--cur_position;
		}
		std::reverse(result[i - q].begin(), result[i - q].end());
	}

	return result;
}

void out(const std::vector<int>& poly) {
	for (size_t it = 0; it != poly.size(); ++it) {
		if (poly[it] != 0) {
			if (!((it != 0) && (poly[it] == 1)))
				std::cout << poly[it];
			if (it != 0)
				std::cout << "x^" << it;
			if (it != (poly.size() - 1))
				std::cout << " + ";
		}
		if ((poly[it] == 0) && (poly.size() == 1))
			std::cout << poly[it];
	}
}

std::vector<int> sum(std::vector<int>& poly1, std::vector<int>& poly2, int n) {
	std::vector<int> max_poly;
	if (poly1.size() > poly2.size())
		max_poly = poly1;
	else
		max_poly = poly2;
	size_t minimum = std::min(poly1.size(), poly2.size());
	size_t maximum = std::max(poly1.size(), poly2.size());
	std::vector<int> result(maximum, 0);
	for (size_t it = 0; it != minimum; ++it) {
		result[it] = (poly1[it] + poly2[it]) % n;
	}
	for (size_t it = minimum; it != maximum; ++it) {
		result[it] = max_poly[it];
	}

	while ((result.size() > 1) && (result[result.size() - 1] == 0))
		result.resize(result.size() - 1);

	return result;
}

std::vector<int> subtraction(std::vector<int>& poly1, std::vector<int>& poly2, int n) {
	std::vector<int> max_poly;
	if (poly1.size() > poly2.size())
		max_poly = poly1;
	else
		max_poly = poly2;
	size_t minimum = std::min(poly1.size(), poly2.size());
	size_t maximum = std::max(poly1.size(), poly2.size());
	std::vector<int> result(maximum, 0);
	for (size_t it = 0; it != minimum; ++it) {
		result[it] = (poly1[it] - poly2[it]) % n;
		if (result[it] < 0)
			result[it] += n;
	}
	for (size_t it = minimum; it != maximum; ++it) {
		result[it] = max_poly[it];
	}

	while (result[result.size() - 1] == 0)
		result.resize(result.size() - 1);
	return result;
}

std::vector<int> multiplication(std::vector<int>& poly1, std::vector<int>& poly2, int n) {
	std::vector<int> result(poly1.size() + poly2.size(), 0);
	std::vector<std::vector<int> > intermediary(poly1.size(), std::vector<int>(poly1.size() + poly2.size(), 0));

	for (size_t i = 0; i != poly1.size(); ++i)
		for (size_t j = 0; j != poly2.size(); ++j)
			intermediary[i][i + j] = (poly1[i] * poly2[j]) % n;

	for (size_t i = 0; i != intermediary.size(); ++i)
		result = sum(result, intermediary[i], n);

	while ((result.size() > 1) && (result[result.size() - 1] == 0))
		result.resize(result.size() - 1);

	return result;
}

std::vector<int> division(std::vector<int>& poly1, std::vector<int>& poly2, int n) {
	std::vector<int> result(poly1.size() - poly2.size() + 1, 0);
	std::vector<int> intermediary = poly1;
	std::vector<int> subtrahend;
	while (intermediary.size() >= poly2.size()) {
		int coeff1 = intermediary[intermediary.size() - 1];
		int coeff2 = poly2[poly2.size() - 1];
		int res_coeff = 0;
		while (((res_coeff*coeff2) % n) != coeff1)
			++res_coeff;
		std::vector<int> current(intermediary.size() - poly2.size() + 1, 0);
		current[intermediary.size() - poly2.size()] = res_coeff;
		result[intermediary.size() - poly2.size()] = res_coeff;
		subtrahend = multiplication(current, poly2, n);
		intermediary = subtraction(intermediary, subtrahend, n);
		while (intermediary[intermediary.size() - 1] == 0)
			intermediary.resize(intermediary.size() - 1);
	}

	return result;
}

std::vector<int> factorization(std::vector<int>& poly1, std::vector<int>& poly2, int n) {
	if (poly1.size() < poly2.size())
		return poly1;
	std::vector<int> result(poly1.size() - poly2.size() + 1, 0);
	std::vector<int> intermediary = poly1;
	std::vector<int> subtrahend;
	if (intermediary.size() < poly2.size())
		return intermediary;
	while (intermediary.size() >= poly2.size()) {
		int coeff1 = intermediary[intermediary.size() - 1];
		int coeff2 = poly2[poly2.size() - 1];
		int res_coeff = 0;
		while (((res_coeff*coeff2) % n) != coeff1)
			++res_coeff;
		std::vector<int> current(intermediary.size() - poly2.size() + 1, 0);
		current[intermediary.size() - poly2.size()] = res_coeff;
		result[intermediary.size() - poly2.size()] = res_coeff;
		subtrahend = multiplication(current, poly2, n);
		intermediary = subtraction(intermediary, subtrahend, n);
		while (intermediary[intermediary.size() - 1] == 0)
			intermediary.resize(intermediary.size() - 1);
	}

	return intermediary;
}


std::vector<int> Eratosthenes(int n) {
	std::vector<bool> prime(n + 1, true);
	prime[0] = prime[1] = false;
	for (int i = 2; i <= n; ++i)
		if (prime[i])
			if (i * 1ll * i <= n)
				for (int j = i*i; j <= n; j += i)
					prime[j] = false;

	std::vector<int> result;
	for (size_t i = 2; i != prime.size(); ++i) {
		if (prime[i])
			result.push_back(i);
	}
	return result;
}

std::vector<int> polypow(std::vector<int>& poly, std::vector<int>& factor, int degree, int mod) {
	std::vector<int> result;
	if (degree == 0) {
		std::vector<int> inter;
		inter.push_back(1);
		result = inter;
		return result;
	}
	if (degree % 2 == 1) {
		result = polypow(poly, factor, degree - 1, mod);
		result = multiplication(result, poly, mod);
		result = factorization(result, factor, mod);
		return result;
	}
	else {
		result = polypow(poly, factor, degree / 2, mod);
		result = multiplication(result, result, mod);
		result = factorization(result, factor, mod);
		return result;
	}
}

int gcd(int a, int b) {
	if (b == 0)
		return a;
	else
		return gcd(b, a % b);
}


std::vector<std::vector<int> > find_primitive(std::vector<std::vector<int> >& field, std::vector<int>& irreducible, int mod) {
	std::vector<std::vector<int> > result;
	std::vector<int> neutral(1, 1);
	int group_size = field.size() - 1;
	std::vector<int> prime = Eratosthenes(group_size);
	for (size_t i = 0; i != group_size; ++i) {
		int divisor_count = 0;
		int count = 0;
		std::vector<int> element = field[i];
		for (size_t j = 0; j != prime.size(); ++j) {
			if ((group_size % prime[j]) == 0) {
				++divisor_count;
				size_t deg = group_size / prime[j];
				std::vector<int> v = polypow(element, irreducible, deg, mod);
				if (v != neutral)
					++count;
			}
		}
		if (count == divisor_count) {
			result.push_back(element);
			break;
		}
	}

	for (size_t index = 2; index <= group_size; ++index) {
		int common_divisor = gcd(index, (group_size));
		if (common_divisor == 1) {
			std::vector<int> polynom = polypow(result[0], irreducible, index, mod);
			result.push_back(polynom);
		}
	}
	
	return result;
}

std::vector<std::vector<int> > find_generator(std::vector<std::vector<int> >& field, std::vector<int>& irreducible, int mod, int degree) {
	std::vector<std::vector<int> > result;
	for (size_t index = 0; index < (field.size() - mod); ++index) {
		std::vector<int> element = field[index];
		int d = 0;
		std::vector<int> polynom = polypow(element, irreducible, mod, mod);
		for (d = 2; d <= degree; ++d) {
			int cur_degree = binpow(mod, d);
			std::vector<int> polynom = polypow(element, irreducible, cur_degree, mod);
			if (polynom == element)
				break;
		}
		if (d == degree)
			result.push_back(element);
	}
	return result;
}


int _tmain(int argc, _TCHAR* argv[])
{
	int degree, mod;
	std::cout << "Degree: ";
	std::cin >> degree;
	std::cout << std::endl << "Mod: ";
	std::cin >> mod;
	std::cout << std::endl;

	std::vector<std::vector<int> > answer = make_polynomials(degree, mod);
	std::set<std::vector<int> > all_poly;
	for (size_t i = 0; i != answer.size(); ++i)
		all_poly.insert(answer[i]);

	std::vector<std::vector<std::vector<int> > > multiplying_poly;
	for (int i = 0; i != degree; ++i) {
		multiplying_poly.push_back(make_polynomials(i, mod));
	}

	std::set <std::vector<int> > reducible_poly;
	int med = degree / 2;
	for (int i = 1; i <= med; ++i) {
	    for (size_t j = 0; j != multiplying_poly[i].size(); ++j) {
			for (size_t k = 0; k != multiplying_poly[degree - i].size(); ++k) {
				std::vector<int> mult_1 = multiplying_poly[i][j];
				std::vector<int> mult_2 = multiplying_poly[degree - i][k];
				reducible_poly.insert(multiplication(mult_1, mult_2, mod));
			}
		}
	}


	for (const auto& poly : reducible_poly) {
		all_poly.erase(poly);
	}


	std::vector<std::vector<int> > irreducible;
	for (auto it = all_poly.begin(); it != all_poly.end(); ++it) {
		if ((*it)[(*it).size() - 1] == 1) {
			//out(*it);
			irreducible.push_back(*it);
			//std::cout << std::endl;
		}
	}

	std::cout << "Irreducible polynomials: " << irreducible.size() << std::endl;
	for (size_t it = 0; it != irreducible.size(); ++it) {
		std::cout << it << ": ";
		out(irreducible[it]);
		std::cout << std::endl;
	}

	size_t number;
	std::cout << "Choose one of these irreducible polynomials: ";
	std::cin >> number;
	while (number >= irreducible.size()) {
		std::cout << std::endl << "This number is incorrect. Write another one: ";
		std::cin >> number;
	}
	out(irreducible[number]);
	std::cout << std::endl;
	
	std::vector<std::vector<int> > field;
	int current_degree = irreducible[number].size() - 2;
	while (current_degree > 0) {
		std::vector<std::vector<int> > remainders = make_polynomials(current_degree, mod);
		for (size_t it = 0; it != remainders.size(); ++it)
			field.push_back(remainders[it]);
		--current_degree;
	}
	int current_mod = mod - 1;
	while (current_mod >= 0) {
		std::vector<int> numbers;
		numbers.push_back(current_mod);
		field.push_back(numbers);
		--current_mod;
	}

	std::cout << std::endl << "Elements of the field: " << field.size() << std::endl;
	for (size_t i = 0; i != field.size(); ++i) {
		std::cout << i << ": ";
		out(field[i]);
		std::cout << std::endl;
	}

	

	int operation;
	size_t group_size = field.size() - 1;
	std::cout << std::endl << "Let's compute sums, products and quotients!";
	std::cout << std::endl << "Sums: enter 1";
	std::cout << std::endl << "Products: enter 2";
	std::cout << std::endl << "Quotients: enter 3";
	std::cout << std::endl << "To continue: enter 0";
	std::cout << std::endl;
	std::cin >> operation;
	while (operation != 0) {
		if (operation == 1) {
			int num1, num2;
			std::cout << "Please, enter the number of the first component: ";
			std::cin >> num1;
			std::cout << std::endl;
			out(field[num1]);
			std::cout << std::endl << "Please, enter the number of the second component: ";
			std::cin >> num2;
			std::cout << std::endl;
			out(field[num2]);
			std::cout << std::endl << "Sum: ";
			std::vector<int> s = sum(field[num1], field[num2], mod);
			out(s);
			std::cout << std::endl;
		}

		if (operation == 2) {
			int num1, num2;
			std::cout << "Please, enter the number of the first multiplicand: ";
			std::cin >> num1;
			std::cout << std::endl;
			out(field[num1]);
			std::cout << std::endl << "Please, enter the number of the second multiplicand: ";
			std::cin >> num2;
			std::cout << std::endl;
			out(field[num2]);
			std::cout << std::endl << "Product: ";
			std::vector<int> product = multiplication(field[num1], field[num2], mod);
			product = factorization(product, irreducible[number], mod);
			out(product);
			std::cout << std::endl;
		}

		if (operation == 3) {
			int num1, num2;
			std::cout << "Please, enter the number of the dividend: ";
			std::cin >> num1;
			std::cout << std::endl;
			out(field[num1]);
			std::cout << std::endl << "Please, enter the number of the divisor: ";
			std::cin >> num2;
			std::cout << std::endl;
			out(field[num2]);
			std::cout << std::endl << "Quotient: ";
			std::vector<int> inverse = polypow(field[num2], irreducible[number], (group_size - 1), mod);
			std::vector<int> product = multiplication(field[num1], inverse, mod);
			product = factorization(product, irreducible[number], mod);
			out(product);
			std::cout << std::endl;
		}
		std::cout << std::endl << "Enter the number of operation: ";
		std::cin >> operation;
	}

	std::vector<std::vector<int> > primitive = find_primitive(field, irreducible[number], mod);
	std::cout << std::endl << "Primitive elements: " << primitive.size() << std::endl;
	for (size_t i = 0; i != primitive.size(); ++i) {
		out(primitive[i]);
		std::cout << std::endl;
	}

	std::vector<std::vector<int> > generating = find_generator(field, irreducible[number], mod, degree);

	std::cout << std::endl << "Generating elements: " << generating.size() << std::endl;
	for (size_t i = 0; i != generating.size(); ++i) {
		out(generating[i]);
		std::cout << std::endl;
	}

	std::cout << std::endl << std::endl << "Let's construct all the isomorphisms between two fields." << std::endl;
	std::cout << "Enter the number of the first irreducible polynomial: ";
	size_t number1, number2;
	std::cin >> number1;
	while (number1 >= irreducible.size()) {
		std::cout << std::endl << "This number is incorrect. Write another one: ";
		std::cin >> number1;
	}
	out(irreducible[number1]);
	std::cout << std::endl;

	std::cout << "Enter the number of the second irreducible polynomial: ";
	std::cin >> number2;
	while (number2 >= irreducible.size()) {
		std::cout << std::endl << "This number is incorrect. Write another one: ";
		std::cin >> number2;
	}
	out(irreducible[number2]);
	std::cout << std::endl;
	std::cout << std::endl;

	generating = find_generator(field, irreducible[number2], mod, degree);
	std::vector<std::vector<int> > all_matchings;
	for (size_t iteration = 0; iteration != generating.size(); ++iteration) {

		std::vector<Isomorphism> matching(field.size());
		std::vector<int> ex(2);
		ex[0] = 0;
		ex[1] = 1;
		matching[0].preimage = ex;
		matching[0].image = generating[iteration];

		std::vector<int> phi = generating[iteration];

		std::vector<int> zero(1, 0);
		std::vector<int> insertion;
		for (size_t it = 0; it != irreducible[number1].size(); ++it) {
			std::vector<int> iter(1, irreducible[number1][it]);
			std::vector<int> result = polypow(phi, irreducible[number2], it, mod);
			result = multiplication(iter, result, mod);
			insertion = sum(insertion, result, mod);
		}


		std::vector<int> remainder = factorization(insertion, irreducible[number2], mod);
		if (remainder == zero) {
			matching[0].image = generating[iteration];
			size_t match_it = 1;

			for (size_t index = 0; index != field.size(); ++index) {
				if (field[index] != ex) {
					matching[match_it].preimage = field[index];
					++match_it;
				}
			}

			for (match_it = 1; match_it != matching.size(); ++match_it) {
				std::vector<int> current_preimage = matching[match_it].preimage;
				std::vector<int> current_image;
				for (size_t i = 0; i != current_preimage.size(); ++i) {
					std::vector<int> iter(1, current_preimage[i]);
					std::vector<int> result = polypow(phi, irreducible[number2], i, mod);
					result = multiplication(iter, result, mod);
					current_image = sum(current_image, result, mod);
				}
				matching[match_it].image = current_image;
			}

			for (size_t i = 0; i != matching.size(); ++i) {
				out(matching[i].preimage);
				std::cout << " -> ";
				out(matching[i].image);
				std::cout << std::endl;
			}
			std::cout << std::endl << std::endl;
		}
	}
	std::cout << std::endl;
	system("pause");
	return 0;
}

