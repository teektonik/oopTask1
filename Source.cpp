#include <iostream>
#include <vector>
#include <stdio.h>
#include <map>
#include <cmath>
#include <numeric>
#include <string>
#include <functional>
#include <ostream>
#include <iomanip>
#include <sstream>

using namespace std;

class BigInt {
public:
	static const int BASE = 1000000000;
	// внутреннее хранилище числа
	std::vector<int> digits;

	// знак числа
	bool is_negative;
	void _remove_leading_zeros();
	BigInt();
	BigInt(int);
	BigInt(std::string); // бросать исключение std::invalid_argument при ошибке
	BigInt(const BigInt&);
	~BigInt() = default;

	BigInt& operator=(const BigInt&);  //возможно присваивание самому себе!

	BigInt operator~() const;

	BigInt& operator++();
	const BigInt operator++(int) const;
	BigInt& operator--();
	const BigInt operator--(int) const;

	BigInt& operator+=(const BigInt&);
	BigInt& operator*=(const BigInt&);
	BigInt& operator-=(const BigInt&);
	BigInt& operator/=(const BigInt&);
	BigInt& operator^=(const BigInt&);
	BigInt& operator%=(const BigInt&);
	BigInt& operator&=(const BigInt&);
	BigInt& operator|=(const BigInt&);

	BigInt operator+() const;  // unary +
	BigInt operator-() const;  // unary -
	void toDecimal(std::vector<bool>&);
	bool operator==(const BigInt&) const;
	bool operator!=(const BigInt&) const;
	bool operator<(const BigInt&) const;
	bool operator>(const BigInt&) const;
	bool operator<=(const BigInt&) const;
	bool operator>=(const BigInt&) const;
	void _shift_right();
	/*operator int() const;
	operator std::string() const;*/
	bool even() const;
	bool odd() const;
	size_t size() const;  // size in bytes
};
size_t BigInt::size() const
{
	return this->digits.size();
}
std::ostream& operator<<(std::ostream& o, const BigInt& i) {
	if (i.digits.empty())
		o << 0;
	else {
		if (i.is_negative)
			o << '-';
		o << i.digits.back();
		char old_fill = o.fill('0');
		for (long long j = static_cast<long long>(i.digits.size()) - 2; j >= 0; --j)
			o << std::setw(9) << i.digits[j];
		o.fill(old_fill);
	}
	return o;
}

void BigInt::_remove_leading_zeros()
{
	while (this->digits.size() > 1 && this->digits.back() == 0) {
		this->digits.pop_back();
	}

	if (this->digits.size() == 1 && this->digits[0] == 0) this->is_negative = false;
}

BigInt::BigInt()
{
	this->is_negative = false;
}

BigInt::BigInt(int a)
{
	if (a < 0) { this->is_negative = true; a = -a; }
	else
		this->is_negative = false;
	do {
		this->digits.push_back(a % BigInt::BASE);
		a /= BigInt::BASE;
	} while (a != 0);
}

BigInt::BigInt(std::string str)
{
	

	if (str.length() == 0) {
		this->is_negative = false;

	}
	else {
		if (str[0] == '-') {
			this->is_negative = true;
			str = str.substr(1);
		}
		else {
			this->is_negative = false;
		}


		for (long long i = str.length(); i > 0; i -= 9) {
			if (i < 9)
				this->digits.push_back(atoi(str.substr(0, i).c_str()));
			else
				this->digits.push_back(atoi(str.substr(i - 9, 9).c_str()));
		}
	}
	this->_remove_leading_zeros();
}

BigInt::BigInt(const BigInt& a)
{
	*this = a;

}

BigInt& BigInt::operator=(const BigInt& a)
{

	if (a == *this)
		return *this;
	this->digits.resize(a.digits.size());
	this->digits = a.digits;
	this->is_negative = a.is_negative;
	return *this;

}

BigInt BigInt::operator~() const
{
	return BigInt();
}


BigInt BigInt::operator+() const
{
	BigInt a = *this;
	a.is_negative = false;
	return a;
	
}

BigInt BigInt::operator-() const
{
	BigInt a = *this;
	a.is_negative = !a.is_negative;
	return a;
}

bool BigInt::operator==(const BigInt& a) const
{
	if (this->is_negative != a.is_negative) {
		return false;
	}
	if (a.size() != this->size())
		return false;
	if ((a.digits.empty() || this->digits.empty()) && a.digits[0] == 0 && this->digits[0] == 0 && a.size() == 1 && this->size() == 1)
		return true;
	for (int i = 0; i < a.size(); i++)
	{
		if ((a.digits[i] != this->digits[i]) == true)
			return false;
	}
	return true;
}

bool BigInt::operator!=(const BigInt& a) const
{
	return !(a == *this);
}

bool BigInt::operator>(const BigInt& a) const
{
	if (a == *this)
		return false;
	if (a.is_negative)
	{
		if (this->is_negative)
			(-a < -(*this));
		else
			return true;
	}
	else if (this->is_negative)
		return false;
	else {
		if (a.size() != this->size())
			return a.size() < this->size();
		else {
			for (int i = 0; i < a.size(); i++)
				if (a.digits[i] != this->digits[i])
					return a.digits[i] < this->digits[i];
		}

	}
}

bool BigInt::operator<(const BigInt& a) const
{
	return !(*this > a);
}

bool BigInt::operator<=(const BigInt& a) const
{
	return (*this < a || *this == a);
}

bool BigInt::operator>=(const BigInt& a) const
{
	return (*this > a || *this == a);
}
BigInt& BigInt::operator+=(const BigInt& a)
{

	if (this->is_negative) {
		if (a.is_negative) {
			*this = -*this;
			BigInt tmp = -a;
			*this += tmp;
			*this = -*this;
			return *this;
		}
		*this = -*this;
		*this -= a;
		*this = -*this;

		return *this;
	}
	else if (a.is_negative) {
		BigInt tmp = a;
		tmp.is_negative = false;
		*this -= tmp;
		return *this;
	}

	int carry = 0;
	for (size_t i = 0; i < min(this->digits.size(), a.digits.size()) || carry != 0; ++i) {
		if (i == this->digits.size()) {
			this->digits.push_back(0);
		}
		this->digits[i] += carry + (i < a.digits.size() ? a.digits[i] : 0);
		carry = this->digits[i] >= BigInt::BASE;
		if (carry != 0) {
			this->digits[i] -= BigInt::BASE;
		}
	}

	return *this;
}

BigInt& BigInt::operator-=(const BigInt& a)
{
	if (a == *this) {
		BigInt x(0);
		*this = x;
		return *this;
	}
	if (a.is_negative) {
		BigInt tmp = -a;
		*this += tmp;
		return *this;
	}
	else if (this->is_negative) {
		BigInt tmp = -*this;
		tmp += a;
		*this = -tmp;
		return *this;
	}
	else if (*this < a) {
		BigInt tmp = -*this;
		*this = a;
		*this += tmp;
		*this = -*this;
		return *this;
	}

	int carry = 0;
	for (size_t i = 0; i < a.digits.size(); ++i) {
		this->digits[i] -= carry + a.digits[i];
		carry = this->digits[i] < 0;
		if (carry != 0) {
			this->digits[i] += BigInt::BASE;
		}
	}
	if (carry != 0) {
		this->digits[this->digits.size()-1] += BigInt::BASE;
	}
	this->_remove_leading_zeros();
	return *this;
}
BigInt operator*(const BigInt& a, const BigInt& b) {
	BigInt x;
	x.digits.resize(a.digits.size() + b.digits.size() + 1);
	for (int i = 0; i < a.digits.size(); i++) {
		int carry = 0;
		for (int j = 0; j < b.digits.size() || carry != 0; j++) {
			long long cur = x.digits[i + j] +
				a.digits[i] * 1LL * (j < b.digits.size() ? b.digits[j] : 0) + carry;
			x.digits[i + j] = static_cast<int>(cur % BigInt::BASE);
			carry = static_cast<int>(cur / BigInt::BASE);
		}
	}
	if (a.is_negative != b.is_negative)
		x.is_negative = true;
	x._remove_leading_zeros();
	return x;
}

BigInt& BigInt::operator*=(const BigInt& a)
{
	return *this = *this * a;
}

BigInt operator-(const BigInt& a, const BigInt& b) {
	BigInt tmp = a;
	tmp -= b;
	return tmp;
}

BigInt operator+(const BigInt& a , const BigInt& b) {
	
	BigInt tmp = a;
	tmp += b;
	return tmp;
}


void BigInt::_shift_right() {
	if (this->digits.size() == 0) {
		this->digits.push_back(0);
		return;
	}
	this->digits.push_back(this->digits[this->digits.size() - 1]);
	for (size_t i = this->digits.size() - 2; i > 0; --i)
		this->digits[i] = this->digits[i - 1];
	this->digits[0] = 0;
}
BigInt operator/(const BigInt& a, const BigInt& b) {
	if (b == 0)
		throw std::invalid_argument("Divide by zero");

	BigInt c = b;
	c.is_negative = false;
	BigInt result, current;
	result.digits.resize(a.digits.size());
	for (long long i = static_cast<long long>(a.digits.size()) - 1; i >= 0; --i) {
		current._shift_right();
		current.digits[0] = a.digits[i];
		current._remove_leading_zeros();
		int x = 0, l = 0, r = BigInt::BASE;
		while (l <= r) {
			int m = (l + r) / 2;
			BigInt t = c * BigInt(m);
			if (t <= current) {
				x = m;
				l = m + 1;
			}
			else r = m - 1;
		}

		result.digits[i] = x;
		current -= b * x;
	}

	result.is_negative = a.is_negative != b.is_negative;
	result._remove_leading_zeros();
	return result;
}


BigInt operator%(const BigInt& a, const BigInt& b) {
	BigInt x = a / b;
	BigInt z = x * b;
	BigInt result = a-z;
	if (result.is_negative)
		result += b;
	return result;
}


int divideByTwo(std::vector<int>& a)
{
	int carry = 0, b = 2;
	for (int i = a.size() - 1; i >= 0; --i)
	{
		long long cur = a[i] + carry * 1ll * 1000000000;
		a[i] = static_cast<int>(cur / b);
		carry = static_cast<int>(cur % b);
	}
	while (a.size() > 1 && a.back() == 0)
	{
		a.pop_back();
	}

	return carry;
}

std::vector<int> toBin(BigInt& a)
{
	int count_bits = 0;
	std::vector<int> _currentNumCopy(a.digits);
	std::vector<int> num_2_32Base;

	num_2_32Base.reserve(a.digits.size() * 2);
	while (_currentNumCopy[0] >= 1)
	{
		int carry = divideByTwo(_currentNumCopy);
		if (count_bits / 32 >= num_2_32Base.size())
		{
			num_2_32Base.push_back(carry);
		}
		else
		{
			num_2_32Base[count_bits / 32] ^= (carry << (count_bits % 32));
		}
		count_bits++;
	}

	return num_2_32Base;
}

bool BigInt::odd() const {
	if (this->digits.size() == 0) return false;
	return this->digits[0] & 1;
}

bool BigInt::even() const {
	return !this->odd();
}
BigInt& BigInt::operator^=(const BigInt& a)
{
	BigInt result("1");

	BigInt x(0);
	if (a == x) {
		return result;
	}
	else if (a < x) {
		*this = x;
		return *this;
	}
	
	BigInt b(*this);
	BigInt z(a);
	while (z != 0) {
		if (z.odd()) result *= b;
		b *= b;
		z /= BigInt(2);
	}
	*this = result;
	return *this;
}

BigInt operator^(const BigInt& a, const BigInt& b) {
	BigInt tmp = a;
	tmp ^= b;
	return tmp;
}
void BigInt::toDecimal(std::vector<bool>& binaryRepresentaion)
{
	this->digits.resize(0);
	BigInt base(1000000000);
	for (int i = binaryRepresentaion.size() - 1; i >= 0; i--)
	{
		if (binaryRepresentaion[i] != false)
		{
			BigInt temp(base ^ BigInt(i));
			*this += temp;
		}
	}
}

BigInt operator&(const BigInt& a, const BigInt& b)
{
	vector<int> abin2_32 = toBin(const_cast<BigInt&>(a));
	vector<int> bin_2_32 = toBin(const_cast<BigInt&>(b));

	int minValueSize = min(abin2_32.size(), bin_2_32.size());

	BigInt outputValue;

	for (int i = 0; i < minValueSize; i++)
	{
		outputValue.digits.push_back((i < bin_2_32.size()) ? abin2_32[i] & bin_2_32[i] : 0);
	}

	outputValue._remove_leading_zeros();

	return outputValue;
}

BigInt& BigInt::operator&=(const BigInt& other)
{
	return *this = *this & other;
}

BigInt operator|(const BigInt& a, const BigInt& b)
{
	BigInt outputValue;
	vector<int> abin2_32 = toBin(const_cast<BigInt&>(a));
	vector<int> bin_2_32 = toBin(const_cast<BigInt&>(b));

	int maxValueSize = max(abin2_32.size(), bin_2_32.size());

	for (int i = 0; i < maxValueSize; i++)
	{
		if (i < abin2_32.size() && i >= bin_2_32.size()) 
			outputValue.digits.push_back(bin_2_32[i] | 0);
		else if (i >= abin2_32.size() && i < bin_2_32.size()) 
			outputValue.digits.push_back(bin_2_32[i] | 0);
		else 
			outputValue.digits.push_back(abin2_32[i] | bin_2_32[i]);
	}

	outputValue._remove_leading_zeros();

	return outputValue;
}

BigInt& BigInt::operator|=(const BigInt& other)
{
	return *this = *this | other;
}


BigInt& BigInt::operator/=(const BigInt& a)
{
	return *this = *this / a;
}


BigInt& BigInt::operator%=(const BigInt& a)
{
	return *this = *this % a;
}

BigInt& BigInt::operator++()
{
	return *this += 1;
}

const BigInt BigInt::operator++(int) const
{
	
	return *this - 1;
}

BigInt& BigInt::operator--()
{
	return (*this -= BigInt(1));
}

const BigInt BigInt::operator--(int) const
{
	return BigInt();
}

int main()
{

	
	BigInt b(5);
	BigInt c(126);

	std::cout << "Comparison:" << endl;
	std::cout << "5 > 126" << endl;
	std::cout << "My answer: ";
	std::cout << (b > c) << endl;
	std::cout << "Correct answer: ";
	std::cout << "0" << endl;
	c = -c;
	std::cout << "5 > -126" << endl;
	std::cout << "My answer: ";
	std::cout << (b > c) << endl;
	std::cout << "Correct answer: ";
	std::cout << "1" << endl;
	b = c;
	b = -b;
	std::cout << "126 == -126" << endl;
	std::cout << "My answer: ";
	std::cout << (b == c) << endl;
	std::cout << "Correct answer: ";
	std::cout << "0" << endl;
	c = -c;
	std::cout << "126 == 126" << endl;
	std::cout << "My answer: ";
	std::cout << (b == c) << endl;
	std::cout << "Correct answer: ";
	std::cout << "1" << endl;
	std::cout << "126 >= 126" << endl;
	std::cout << "My answer: ";
	std::cout << (b >= c) << endl;
	std::cout << "Correct answer: ";
	std::cout << "1" << endl;
	BigInt x("1279382781331273788732398745768978");
	BigInt a("12300078658673927675889379123798");
	std::cout << "1279382781331273788732398745768978 < 12300078658673927675889379123798" << endl;
	std::cout << "My answer: ";
	std::cout << (x < a) << endl;
	std::cout << "Correct answer: ";
	std::cout << "0" << endl;
	std::cout << endl;

	std::cout << "Plus:" << endl;
	std::cout << "1279382781331273788732398745768978 + 12300078658673927675889379123798" << endl;
	std::cout << "My answer: ";
	std::cout << (x + a) << endl;
	std::cout << "Correct answer: ";
	std::cout << "1291682859989947716408288124892776" << endl;
	b =  -5;
	std::cout << "126 + (-5)" << endl;
	std::cout << "My answer: ";
	std::cout << (c + b) << endl;
	std::cout << "Correct answer: ";
	std::cout << "121" << endl;
	std::cout << "(-5) + 126" << endl;
	std::cout << "My answer: ";
	std::cout << (b + c) << endl;
	std::cout << "Correct answer: ";
	std::cout << "121" << endl;
	c = -c;
	std::cout << "(-5) + (-126)" << endl;
	std::cout << "My answer: ";
	std::cout << (b + c) << endl;
	std::cout << "Correct answer: ";
	std::cout << "-131" << endl;
	std::cout << endl;

	std::cout << "Minus:" << endl;
	std::cout << "1279382781331273788732398745768978 - 12300078658673927675889379123798" << endl;
	std::cout << "My answer: ";
	std::cout << (x - a) << endl;
	std::cout << "Correct answer: ";
	std::cout << "1267082702672599861056509366645180" << endl;
	c = -c;
	std::cout << "126 - (-5)" << endl;
	std::cout << "My answer: ";
	std::cout << (c - b) << endl;
	std::cout << "Correct answer: ";
	std::cout << "131" << endl;
	std::cout << "(-5) - 126" << endl;
	std::cout << "My answer: ";
	std::cout << (b - c) << endl;
	std::cout << "Correct answer: ";
	std::cout << "-131" << endl;
	c = -c;
	std::cout << "(-5) - (-126)" << endl;
	std::cout << "My answer: ";
	std::cout << (b - c) << endl;
	std::cout << "Correct answer: ";
	std::cout << "121" << endl;
	std::cout << endl;

	std::cout << "Multiply:" << endl;
	std::cout << "1279382781331273788732398745768978 * 12300078658673927675889379123798" << endl;
	std::cout << "My answer: ";
	std::cout << (x * a) << endl;
	std::cout << "Correct answer: ";
	std::cout << "15736508844927693021137653772270497294378591992272917813369938444" << endl;
	c = -c;
	std::cout << "126 * (-5)" << endl;
	std::cout << "My answer: ";
	std::cout << (c * b) << endl;
	std::cout << "Correct answer: ";
	std::cout << "-630" << endl;
	std::cout << "(-5) * 126" << endl;
	std::cout << "My answer: ";
	std::cout << (b * c) << endl;
	std::cout << "Correct answer: ";
	std::cout << "-630" << endl;
	c = -c;
	std::cout << "(-5) * (-126)" << endl;
	std::cout << "My answer: ";
	std::cout << (b * c) << endl;
	std::cout << "Correct answer: ";
	std::cout << "630" << endl;
	std::cout << endl;
	BigInt one("12312312312312");
	BigInt two = 12;

	std::cout << "Divide:" << endl;
	std::cout << "12312312312312 / 12" << endl;
	std::cout << "My answer: ";
	std::cout << (one / two) << endl;
	std::cout << "Correct answer: ";
	std::cout << "1026026026026" << endl;
	c = -c;
	std::cout << "126 / (-5)" << endl;
	std::cout << "My answer: ";
	std::cout << (c / b) << endl;
	std::cout << "Correct answer: ";
	std::cout << "-25" << endl;
	std::cout << "(-5) / 126" << endl;
	std::cout << "My answer: ";
	std::cout << (b / c) << endl;
	std::cout << "Correct answer: ";
	std::cout << "0" << endl;
	c = -c;
	std::cout << "(-126) / (-5)" << endl;
	std::cout << "My answer: ";
	std::cout << (c / b) << endl;
	std::cout << "Correct answer: ";
	std::cout << "25" << endl;
	std::cout << endl;
	BigInt p("132");
	BigInt o("23");
	std::cout << "132 / 23" << endl;
	std::cout << "My answer: ";
	std::cout << (p % o) << endl;
	std::cout << "Correct answer: ";
	std::cout << "17" << endl;
	std::cout << endl;

	std::cout << "Bit operations:" << endl;
	std::cout << "12312312312312 & 12" << endl;
	std::cout << "My answer: ";
	std::cout << (one & two) << endl;
	std::cout << "Correct answer: ";
	std::cout << "8" << endl;
	c = -c;
	one = 123123123;
	std::cout << "123123123 | 12" << endl;
	std::cout << "My answer: ";
	std::cout << (one | two) << endl;
	std::cout << "Correct answer: ";
	std::cout << "123123135" << endl;
	
	std::cout << endl;
	BigInt k("123");
	BigInt l("12");
	std::cout << "Pow:" << endl;
	std::cout << "123 ^ 12" << endl;
	std::cout << "My answer: ";
	std::cout << (k ^ l) << endl;
	std::cout << "Correct answer: ";
	std::cout << "11991163848716906297072721" << endl;
	l = -l;
	std::cout << "123 ^ (-12)" << endl;
	std::cout << "My answer: ";
	std::cout << (k ^ l) << endl;
	std::cout << "Correct answer: ";
	std::cout << "0" << endl;
	
	return 0;
}
