{\rtf1\ansi\ansicpg1252\cocoartf1671\cocoasubrtf600
{\fonttbl\f0\fnil\fcharset0 Menlo-Regular;}
{\colortbl;\red255\green255\blue255;\red36\green38\blue41;\red235\green236\blue237;\red114\green121\blue129;
\red13\green0\blue129;\red43\green39\blue19;}
{\*\expandedcolortbl;;\cssrgb\c18824\c20000\c21176;\cssrgb\c93725\c94118\c94510;\cssrgb\c52157\c54902\c57647;
\cssrgb\c6275\c6275\c58039;\cssrgb\c22353\c20000\c9412;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\sl300\partightenfactor0

\f0\fs26 \cf2 \cb3 \expnd0\expndtw0\kerning0
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]\
\
\pard\pardeftab720\sl300\partightenfactor0
\cf4 #include <RcppArmadillo.h>\cf2 \
\cf4 #include <RcppEigen.h>\cf2 \
\
// [[Rcpp::export]]\
SEXP armaMatMult(arma::mat A, arma::mat B)\{\
    arma::mat C = A * B;\
\
    \cf5 return\cf2  Rcpp::wrap(C);\
\}\
\
// [[Rcpp::export]]\
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B)\{\
    Eigen::MatrixXd C = A * B;\
\
    \cf5 return\cf2  Rcpp::wrap(C);\
\}\
\
// [[Rcpp::export]]\
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B)\{\
    Eigen::MatrixXd C = A * B;\
\
    \cf5 return\cf2  Rcpp::wrap(C);\
\}\cf6 \
}