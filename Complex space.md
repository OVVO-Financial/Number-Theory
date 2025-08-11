# Factors in the complex space

Below is the complex space associated with `n=309`.  Each triangle represents a complex number whose lower [corresponding real number](https://github.com/OVVO-Financial/Number-Theory/blob/master/Number%20Theory%20Papers/i.pdf](https://github.com/OVVO-Financial/Number-Theory/blob/master/Number%20Theory%20Papers/Real_Pair_Representation_of_Complex_Numbers_and_Quaternions.pdf) lies within the interval `[3,sqrt(n)]`.

* Blue triangles are complex numbers whose corresponding real product is less than `n`.

* Red triangles are complex numbers whose corresponding real product is greater than `n`.

* The green triangle represents the complex number associated with the factors of `n=309`, `53+50i` which has corresponding real numbers `[3,103]`.

![Complex Space](https://github.com/OVVO-Financial/Number-Theory/blob/master/Images/Complex%20plane.png)

## Possible complex factors
All possible factors reside along the intersection of red & blue triangles.  We can highlight this strip in green.  The question becomes, ***how do we best navigate this strip?***

If we could define this green strip, we could search only the integer complex numbers along the strip.  Objectively, it looks like a rotated logarithmic function, but, along with every other [prime function appearing logarithmic](https://github.com/OVVO-Financial/Number-Theory/blob/master/Number%20Theory%20Papers/On%20the%20Distribution%20of%20Prime%20Numbers.pdf), it defies fitting!

![Factor Strip](https://github.com/OVVO-Financial/Number-Theory/blob/master/Images/Factor%20Strip%20in%20Green.jpeg)

### Trial Division
Trial division diagonally checks all complex points as illuastrated below, where all of the highlighted complex numbers have a lower real number = 3.  We can see its effectiveness for the upper part of the green strip which becomes closer to 45 degrees.

![Trial Div](https://github.com/OVVO-Financial/Number-Theory/blob/master/Images/Trial%20Division%20by%203.jpeg)

### Fermat
Fermat's difference of squares `N = a^2 - b^2` vertically checks all complex points illustrated below, where all of the highlighted complex numbers are equal to `a` in `b^2 = a^2 - N`.  Fermat is most effective for the lower part of the green strip which is more vertical at its base near `sqrt(n)`.

Fermat's range extends horizontally such that all values of `b` are considered below the corresponding horizontal line (hence the squares!).

#### Vertical Fermat
![Fermat](https://github.com/OVVO-Financial/Number-Theory/blob/master/Images/Vertical%20Fermat.jpeg)

#### Horizontal Fermat
![H Fermat](https://github.com/OVVO-Financial/Number-Theory/blob/master/Images/Horizontal%20Fermat.jpeg)


### Combined methods
The above images explain why combining trial division and Fermat's method is more effective than either on its own (while each has the ability to be sieved, namely primes only for trial division and the [Fermat sieve here](https://github.com/OVVO-Financial/Number-Theory/blob/master/Number%20Theory%20Papers/Fermat%20Sieve%20Using%20Complex%20Numbers.pdf)).  However, we are still left with the middle section of the strip to navigate.  [Complex trial multiplication](https://github.com/OVVO-Financial/Number-Theory/blob/Prime-Factorization/Complex%20Trial%20Multiplication.md) is most effective in this area, permitting us to [combine all 3 methods](https://github.com/OVVO-Financial/Number-Theory/blob/Prime-Factorization/julia/Simultaneous_Complex_Factorization.jl) to simultaneously navigate their most effective part of the factor strip.

The yellow triangles represent the number of steps to completely scan the entire complex space using Fermat & trial division methods simultaneously.

![Simultaneous](https://github.com/OVVO-Financial/Number-Theory/blob/master/Images/Complex%20Space%20Factorization%201.jpeg)

### One final insight
We have one trick to turn this rotated logarithmic factor strip into a straight vertical line in the complex space.  *If we square the complex space*, we have the following image.  The complex factors, when squared, have a `real` part equal to `n`, our number we wish to factor. The `imaginary` coefficient is equal to `2ab`, from our Fermat terms in `a^2 + b^2 = N`.

Thus we only need to test if the square root of these squared complex numbers have integer `real` and `imaginary` coefficients.

#### Original complex space
Our original complex factor for `n=309` is `53+50i`.
![Factor Strip](https://github.com/OVVO-Financial/Number-Theory/blob/master/Images/Factor%20Strip%20in%20Green.jpeg)

#### Squared complex space
The final form of our squared complex factor will be `n+2abi`.  Our squared complex factor for `n=309` is `(53+50i)^2 = 309+5300i`.  The large `imaginary` coefficient may seem daunting, but remember, we only need to test (via `sqrt(n + 2abi)`)<sup>[1](#footnote1)</sup> for integer solutions using intervals of `2ab` which grow quickly.  In our example where:
```
n=309
```
our solutions were:
```
a=53 and b=50 
```
thus
```
2ab=5300
```

Also, the [Fermat sieves](https://github.com/OVVO-Financial/Number-Theory/blob/master/Number%20Theory%20Papers/Fermat%20Sieve%20Using%20Complex%20Numbers.pdf) for `a` and `b` are still applicable.

The green vertical line, formerly our factor strip, is now simply our real number `n=309`!

![Squared complex space](https://github.com/OVVO-Financial/Number-Theory/blob/master/Images/Complex%20squared.png)

All of this yields the final complex polynomial<sup>[2](#footnote2)</sup> to solve (perhaps using a complex Newton-Raphson method):
```
(a + bi)^2 -n -2abi = 0
```
And here is the equivalent polynomial using only real numbers:
```
sqrt(n^2 + [2ab]^2) - (a^2 + b^2) = 0
```

where our intervals for each component are defined by:
```
a interval: 
a = [minimum real, maximimum real]
a = [ceil(sqrt(n)),3+((n-9)/6)]

b interval:
b = [minimum imaginary, maximum imaginary]
b = [1,(n-9)/6] 
```

##### Footnotes:
<a name="footnote1">1</a>: This can be reduced to an integer hypoteneuse test due to the right triangle each complex number represents.  The hypoteneuse of a right triangle with sides `309,5300` is equal to `sqrt(309^2 + 5300^2) = 5309` which is equal to the sum of our squared complex factor components `(50^2 + 53^2)`.  More generally, we can state `sqrt(n^2 + [2ab]^2) = (a^2 + b^2)` where integer solutions for `a and b` will yield the factors of `n`.

<a name="footnote2">2</a>: This is of course equivalent to Fermat `a^2 - b^2 = n`, but the hope is alternative representations yield new insights & methods.

The above plots were generated with the [complex space generator](https://github.com/OVVO-Financial/Number-Theory/blob/master/R/Complex%20Space%20Generator.R) in R.
