
float **PNm(l, theta)
int l;
float theta;
{
	static int prvl = -1, l3;
	static float **P, prvtheta = 3.2;
	int N, m;
	float **farray2();

	if (l == prvl && theta == prvtheta)
		return (P);

	if (prvl != -1)
		free_farray2(P, -1, 1, -l3);
	prvl = l;

	l3 = l;
	if (l3 < 3)
		l3 = 3;
	P = farray2(-1, 1, -l3, l3);
	for (N = -1; N <= 1; N++)
		for (m = -l3; m <= l3; m++)
			P[N][m] = 0.0;

	if (l)
		rotmatrix(1, l, theta, P);
	else
		P[0][0] = 1.0;

	return (P);

}
