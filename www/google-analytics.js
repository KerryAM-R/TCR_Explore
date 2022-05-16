{\rtf1\ansi\ansicpg1252\cocoartf2580
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red51\green110\blue109;\red246\green246\blue246;\red38\green38\blue38;
\red15\green112\blue1;\red169\green14\blue26;\red83\green83\blue83;}
{\*\expandedcolortbl;;\cssrgb\c25098\c50196\c50196;\cssrgb\c97255\c97255\c97255;\cssrgb\c20000\c20000\c20000;
\cssrgb\c0\c50196\c0;\cssrgb\c72941\c12941\c12941;\cssrgb\c40000\c40000\c40000;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs26 \cf2 \cb3 \expnd0\expndtw0\kerning0
// Initial Tracking Code\cf4 \
(\cf5 function\cf4 (i,s,o,g,r,a,m)\{\
  i['\cf6 GoogleAnalyticsObject\cf4 ']=r;\
  i[r]=i[r] || \
  \cf5 function\cf4 ()\{\
    (i[r].q=i[r].q||[]).push(arguments);\
  \},i[r].l=\cf7 1\cf4 *\cf5 new\cf4  Date();\
  a=s.createElement(o),\
  m=s.getElementsByTagName(o)[\cf7 0\cf4 ];\
  a.\cf5 async\cf4 =\cf7 1\cf4 ;\
  a.src=g;\
  m.parentNode.insertBefore(a,m);\
\})(window,document,'\cf6 script\cf4 ',\
  '\cf6 https://www.google-analytics.com/analytics.js\cf4 ','\cf6 ga\cf4 ');\
\
ga('\cf6 create\cf4 ', '\cf6 <YOUR_APP_ID>\cf4 ', '\cf6 auto\cf4 ');\
ga('\cf6 send\cf4 ', '\cf6 pageview\cf4 ');\
\
\cf2 // Event Tracking Code\cf4 \
$(document).on('\cf6 shiny:inputchanged\cf4 ', \cf5 function\cf4 (event) \{\
  \cf5 if\cf4 (event.name == '\cf6 bins\cf4 ' || event.name == '\cf6 col\cf4 ')\{\
    ga('\cf6 send\cf4 ', '\cf6 event\cf4 ', '\cf6 input\cf4 ',\
      '\cf6 updates\cf4 ', event.name, event.value);\
  \}\
\});\
\
\cf2 // User Tracking Code\cf4 \
$(document).one('\cf6 shiny:idle\cf4 ', \cf5 function\cf4 () \{\
  ga('\cf6 set\cf4 ','\cf6 userId\cf4 ', Shiny.user);\
\});}