
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int id_x12,id_xr,id_xs,N,i,j,isC;
    float *res,*x0,*x,mu_x1,mu_x2,mu_xr,mu_xs;
    mxArray *var,*X,*AX,*AmX,*AtmX,*lhs[1],*rhs[4];
    mxArray *Wx1,*Wtx1,*wp1,*Wx2,*Wtx2,*wp2,*Wxr,*Wtxr,*wpr,*Wxs,*Wtxs,*wps;

    if(mxIsComplex(prhs[0]))
    {isC=1;}
    else
    {isC=0;}
    id_x12=(int)mxGetScalar(mxGetField(prhs[1],0,"id_x12"));
    if(id_x12==1)
    {   N=(int)mxGetNumberOfElements(prhs[0])/2;
        mu_x1=(float)mxGetScalar(mxGetField(prhs[1],0,"mu_x1"));
        mu_x2=(float)mxGetScalar(mxGetField(prhs[1],0,"mu_x2"));
        if(isC==1)
        {   plhs[0]=mxCreateNumericMatrix(2*N,1,mxSINGLE_CLASS,1);
            X=mxCreateNumericMatrix(N,1,mxSINGLE_CLASS,1);
        }
        else
        {   plhs[0]=mxCreateNumericMatrix(2*N,1,mxSINGLE_CLASS,0);
            X=mxCreateNumericMatrix(N,1,mxSINGLE_CLASS,0);
        }
        Wx1=(mxArray *)mxGetField(prhs[1],0,"Wx1");
        Wtx1=(mxArray *)mxGetField(prhs[1],0,"Wtx1");
        wp1=(mxArray *)mxGetField(prhs[1],0,"wp1");
        Wx2=(mxArray *)mxGetField(prhs[1],0,"Wx2");
        Wtx2=(mxArray *)mxGetField(prhs[1],0,"Wtx2");
        wp2=(mxArray *)mxGetField(prhs[1],0,"wp2");
    }
    else
    {   N=(int)mxGetNumberOfElements(prhs[0]);
        if(isC==1)
        {plhs[0]=mxCreateNumericMatrix(N,1,mxSINGLE_CLASS,1);}
        else
        {plhs[0]=mxCreateNumericMatrix(N,1,mxSINGLE_CLASS,0);}
    }
    if(mxIsComplex(prhs[0]))
    {AX=mxCreateNumericMatrix(N,1,mxSINGLE_CLASS,1);}
    else
    {AX=mxCreateNumericMatrix(N,1,mxSINGLE_CLASS,0);}

    id_xr=(int)mxGetScalar(mxGetField(prhs[1],0,"id_xr"));
    if(id_xr==1)
    {   mu_xr=(float)mxGetScalar(mxGetField(prhs[1],0,"mu_xr"));
        Wxr=(mxArray *)mxGetField(prhs[1],0,"Wxr");
        Wtxr=(mxArray *)mxGetField(prhs[1],0,"Wtxr");
        wpr=(mxArray *)mxGetField(prhs[1],0,"wpr");
    }
    id_xs=(int)mxGetScalar(mxGetField(prhs[1],0,"id_xs"));
    if(id_xs==1)
    {   mu_xs=(float)mxGetScalar(mxGetField(prhs[1],0,"mu_xs"));
        Wxs=(mxArray *)mxGetField(prhs[1],0,"Wxs");
        Wtxs=(mxArray *)mxGetField(prhs[1],0,"Wtxs");
        wps=(mxArray *)mxGetField(prhs[1],0,"wps");
    }

    AmX=(mxArray *)mxGetField(prhs[1],0,"AmX");
    AtmX=(mxArray *)mxGetField(prhs[1],0,"AtmX");
    var=(mxArray *)mxGetField(prhs[1],0,"var_AtA");

    if(id_x12==1)
    {   x0=mxGetData(prhs[0]);
        x=mxGetData(X);
        for(i=0,j=N;i<N;i++,j++)
        {x[i]=x0[i]+x0[j];}
        if(isC==1)
        {   x0=mxGetImagData(prhs[0]);
            x=mxGetImagData(X);
            for(i=0,j=N;i<N;i++,j++)
            {x[i]=x0[i]+x0[j];}
        }
        rhs[0]=X;lhs[0]=AX;
    }
    else
    {rhs[0]=(mxArray *)prhs[0];lhs[0]=AX;}
    rhs[1]=AmX;rhs[2]=AtmX;rhs[3]=var;
    mexCallMATLAB(1,lhs,4,rhs,"AtAmX");

    res=mxGetData(plhs[0]);
    x=mxGetData(lhs[0]);
    for(i=0;i<N;i++)
    {res[i]=x[i];}
    if(isC==1)
    {   res=mxGetImagData(plhs[0]);
        x=mxGetImagData(lhs[0]);
        for(i=0;i<N;i++)
        {res[i]=x[i];}
    }

    if(id_x12==1)
    {   //update x1+x2
        if(id_xs==1)
        {   lhs[0]=AX;rhs[0]=X;rhs[1]=Wxs;rhs[2]=Wtxs;rhs[3]=wps;
            mexCallMATLAB(1,lhs,4,rhs,"AtAmX");
            res=mxGetData(plhs[0]);
            x=mxGetData(lhs[0]);
            for(i=0;i<N;i++)
            {res[i]+=mu_xs*x[i];}
            if(isC==1)
            {   res=mxGetImagData(plhs[0]);
                x=mxGetImagData(lhs[0]);
                for(i=0;i<N;i++)
                {res[i]+=mu_xs*x[i];}
            }
        }
        if(id_xr==1)
        {   lhs[0]=AX;rhs[0]=X;rhs[1]=Wxr;rhs[2]=Wtxr;rhs[3]=wpr;
            mexCallMATLAB(1,lhs,4,rhs,"AtAmX");
            res=mxGetData(plhs[0]);
            x=mxGetData(lhs[0]);
            for(i=0;i<N;i++)
            {res[i]+=mu_xr*x[i];}
            if(isC==1)
            {   res=mxGetImagData(plhs[0]);
                x=mxGetImagData(lhs[0]);
                for(i=0;i<N;i++)
                {res[i]+=mu_xr*x[i];}
            }
        }
        res=mxGetData(plhs[0]);
        for(i=0,j=N;i<N;i++,j++)
        {res[j]=res[i];}
        if(isC==1)
        {   res=mxGetImagData(plhs[0]);
            for(i=0,j=N;i<N;i++,j++)
            {res[j]=res[i];}
        }
        //update x1
        x0=mxGetData(prhs[0]);
        x=mxGetData(X);
        for(i=0;i<N;i++)
        {x[i]=x0[i];}
        if(isC==1)
        {   x0=mxGetImagData(prhs[0]);
            x=mxGetImagData(X);
            for(i=0;i<N;i++)
            {x[i]=x0[i];}
        }

        lhs[0]=AX;rhs[0]=X;rhs[1]=Wx1;rhs[2]=Wtx1;rhs[3]=wp1;
        mexCallMATLAB(1,lhs,4,rhs,"AtAmX");

        res=mxGetData(plhs[0]);
        x=mxGetData(lhs[0]);
        for(i=0;i<N;i++)
        {res[i]+=mu_x1*x[i];}
        if(isC==1)
        {   res=mxGetImagData(plhs[0]);
            x=mxGetImagData(lhs[0]);
            for(i=0;i<N;i++)
            {res[i]+=mu_x1*x[i];}
        }
        //update x2
        x0=mxGetData(prhs[0]);
        x=mxGetData(X);
        for(i=0,j=N;i<N;i++,j++)
        {x[i]=x0[j];}
        if(isC==1)
        {   x0=mxGetImagData(prhs[0]);
            x=mxGetImagData(X);
            for(i=0,j=N;i<N;i++,j++)
            {x[i]=x0[j];}
        }

        lhs[0]=AX;rhs[0]=X;rhs[1]=Wx2;rhs[2]=Wtx2;rhs[3]=wp2;
        mexCallMATLAB(1,lhs,4,rhs,"AtAmX");

        res=mxGetData(plhs[0]);
        x=mxGetData(lhs[0]);
        for(i=0,j=N;i<N;i++,j++)
        {res[j]+=mu_x2*x[i];}
        if(isC==1)
        {   res=mxGetImagData(plhs[0]);
            x=mxGetImagData(lhs[0]);
            for(i=0,j=N;i<N;i++,j++)
            {res[j]+=mu_x2*x[i];}
        }
        mxDestroyArray(X);
    }
    else
    {   //update x
        if(id_xs==1)
        {   lhs[0]=AX;rhs[0]=(mxArray *)prhs[0];rhs[1]=Wxs;rhs[2]=Wtxs;rhs[3]=wps;
            mexCallMATLAB(1,lhs,4,rhs,"AtAmX");
            res=mxGetData(plhs[0]);
            x=mxGetData(lhs[0]);
            for(i=0;i<N;i++)
            {res[i]+=mu_xs*x[i];}
            if(isC==1)
            {   res=mxGetImagData(plhs[0]);
                x=mxGetImagData(lhs[0]);
                for(i=0;i<N;i++)
                {res[i]+=mu_xs*x[i];}
            }
        }
        if(id_xr==1)
        {
            lhs[0]=AX;rhs[0]=(mxArray *)prhs[0];rhs[1]=Wxr;rhs[2]=Wtxr;rhs[3]=wpr;
            mexCallMATLAB(1,lhs,4,rhs,"AtAmX");
            res=mxGetData(plhs[0]);
            x=mxGetData(lhs[0]);
            for(i=0;i<N;i++)
            {res[i]+=mu_xr*x[i];}
            if(isC==1)
            {   res=mxGetImagData(plhs[0]);
                x=mxGetImagData(lhs[0]);
                for(i=0;i<N;i++)
                {res[i]+=mu_xr*x[i];}
            }
        }
    }
    mxDestroyArray(AX);
}
