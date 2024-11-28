/* given an positive integer n return min((2^N) >= n) */
#include <stdio.h>
#include <stdlib.h>
#include "messageiolib.h"

int min2n(n)
int n;
{
int i;
	if(n<=0) STOP("min2n: illegal argument");
	
	i=2;
	while(i<n) i*=2;
	return i;
}
