    // Transpose gsource_pass to comply with regressit
    transpose(t_gsource_pass,gsource_pass,rsize,n);

    double *t_gsource_pass_fm;
    ZALLOC(t_gsource_pass_fm, rsize_pass*n, double);
    int fm, fma;
    for(fm = 0; fm < n; fm++){
     for(fma = 0; fma < rsize_pass; fma++){
      t_gsource_pass_fm[fm*rsize_pass+fma] = t_gsource_pass[fm*rsize+fma];
     }
    }

    double *gsource_pass_fm;
    ZALLOC(gsource_pass_fm, n*rsize_pass, double);
    for(fm = 0; fm < rsize_pass; fm++){
     for(fma = 0; fma < n; fma++){
      gsource_pass_fm[fm*n+fma] = gsource_pass[fm*n+fma];
     }
    }

    regressit(regans, t_gsource_pass_fm, gtarget, n, rsize_pass) ; //run regression
    mulmat(www, regans, gsource_pass_fm,  1, rsize_pass, n) ; //multiply regans and gsource_pass
    free(t_gsource_pass_fm);
    free(gsource_pass_fm);
//    regressit(regans, t_gsource_pass, gtarget, n, rsize_pass) ; //run regression
//    mulmat(www, regans, gsource_pass,  1, rsize_pass, n) ; //multiply regans and gsource_pass
    /* End of fix */


