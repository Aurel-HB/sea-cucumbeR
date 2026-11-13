// Copyright (C) 2013-2015 Kasper Kristensen
 // License: GPL-2
 
 namespace autodiff {
   /* Gradient */
   template<class Functor, class Type>
   struct gradient_t
   {
     Functor f;
     gradient_t(Functor f_) : f(f_) {}
     AD<Type> userfun(vector<AD<Type> > x){
       return f(x);
     }
     vector<Type> operator()(vector<Type> x0){
       CppAD::vector<AD<Type> > x( x0 );
       CppAD::vector<AD<Type> > y( 1 );
       CppAD::Independent(x);
       y[0] = userfun(x);
       CppAD::ADFun<Type> F(x, y);
       CppAD::vector<Type> x_eval(x0);
       return F.Jacobian(x_eval);
     }
   };
 
   template<class Functor, class Type>
   vector<Type> gradient(Functor F, vector<Type> x){
 #ifdef CPPAD_FRAMEWORK
     gradient_t<Functor, Type> f(F);
     return f(x);
 #endif
 #ifdef TMBAD_FRAMEWORK
     TMBad::ADFun<> G(TMBad::StdWrap<Functor,vector<TMBad::ad_aug> >(F), x);
     G = G.JacFun();
     return G(x);
 #endif
   }
 
   /* Hessian */
   template<class Functor, class Type>
   struct hessian_t
   {
     Functor f;
     gradient_t<Functor,AD<Type> > gr;
     hessian_t(Functor f_) : f(f_), gr(f_) {}
     vector<AD<Type> > userfun(vector<AD<Type> > x){
       return gr(x);
     }
     matrix<Type> operator()(vector<Type> x0){
       CppAD::vector<AD<Type> > x( x0 );
       CppAD::vector<AD<Type> > y( x0 );
       CppAD::Independent(x);
       y = userfun(x);
       CppAD::ADFun<Type> F(x, y);
       CppAD::vector<Type> x_eval(x0);
       vector<Type> ans = F.Jacobian(x_eval);
       return asMatrix(ans, x.size(), x.size());
     }
   };
 
   template<class Functor, class Type>
   matrix<Type> hessian(Functor f, vector<Type> x){
 #ifdef CPPAD_FRAMEWORK
     hessian_t<Functor, Type> H(f);
     return H(x);
 #endif
 #ifdef TMBAD_FRAMEWORK
     TMBad::ADFun<> F(TMBad::StdWrap<Functor,vector<TMBad::ad_aug> >(f), x);
     F = F.JacFun().JacFun();
     vector<Type> Fx = F(x);
     return asMatrix(Fx, x.size(), x.size());
 #endif
   }
 
   /* Jacobian */
   template<class Functor, class Type>
   struct jacobian_t
   {
     Functor f;
     jacobian_t(Functor f_) : f(f_) {}
     vector<AD<Type> > userfun(vector<AD<Type> > x){
       return f(x);
     }
     matrix<Type> operator()(vector<Type> x0){
       CppAD::vector<AD<Type> > x( x0 );
       CppAD::Independent(x);
       CppAD::vector<AD<Type> > y = userfun(x);
       CppAD::ADFun<Type> F(x, y);
       CppAD::vector<Type> x_eval(x0);
       vector<Type> ans = F.Jacobian(x_eval); // By row !
       return asMatrix(ans, x.size(), y.size()).transpose();
     }
   };
 
   template<class Functor, class Type>
   matrix<Type> jacobian(Functor f, vector<Type> x){
 #ifdef CPPAD_FRAMEWORK
     jacobian_t<Functor, Type> J(f);
     return J(x);
 #endif
 #ifdef TMBAD_FRAMEWORK
     TMBad::ADFun<> J(TMBad::StdWrap<Functor,vector<TMBad::ad_aug> >(f), x);
     int m = J.Range();
     J = J.JacFun();
     vector<Type> Jx = J(x);
     return asMatrix(Jx, x.size(), m).transpose();
 #endif
   }
 }