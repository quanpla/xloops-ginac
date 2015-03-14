# R Function #
From Khiem P. H. thesis, we found one implement of R function as following (C++ code).

```
ex fn_R_m1e(const ex &epsilon, const ex &x, const ex &y, const ex &z){
	/*
		Function R_(-1-e) [1, 1, eps; x, y, z]
	*/
	ex R_m1e = 0;
	ex logx = log(x), logy = log(y);

	/* line 1 */
	R_m1e = (logx-logy)/(x-y)	+	epsilon*(logx-logy)/(x-y);

	/* orther lines */
	R_m1e += ( epsilon/(x-y) )/
			(
				pow(logy, 2)/2.0 - pow(logx, 2)/2.0 + Li2(1.0 - z/y) - Li2(1.0 - z/x)
				+ logy*(  fn_eta(z-y, 1.0/(1.0-y) ) - fn_eta(z-y, -1.0/y)  )
				- logx*(  fn_eta(z-x, 1.0/(1.0-x) ) - fn_eta(z-x, -1.0/x)  )
				+ log(1.0-z/y)*fn_eta(z, 1.0/y)  -  log(1.0-z/x)*fn_eta(z, 1.0/x)
			);

	return R_m1e;
}
```