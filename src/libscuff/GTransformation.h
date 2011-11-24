/*
 * GTransformation.h  -- a very simple data structure for 
 *                    -- handling geometric transformations
 *                        
 * homer reid         -- 11/2011
 */

#ifndef GTRANSFORMATION_H
#define GTRANSFORMATION_H

#define GTRANSFORMATION_DISPLACEMENT 0
#define GTRANSFORMATION_ROTATION     1

typedef struct GTransformation
 { 
   int Type;
   double DX[3];
   double M[3][3];
 } GTransformation;

#endif // #ifndef GTRANSFORMATION_H
