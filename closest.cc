#include "closest.h"
#include "point.h"
#include "vector.h"


point closest_point_triangle(const point &A, const point& B, const point& C)
{  

  vector AB = B-A;
  vector AC = C-A;
  vector AP = -A;
  float D1 = dot(AB,AP);
  float D2 = dot(AC,AP);
  vector BP = -B;
  float D3 = dot(AB,BP);
  float D4 = dot(AC,BP);
  float v = D1/(D1-D3);
  vector CP = -C; 
  float D5 = dot(AB,CP);
  float D6 = dot(AC,CP);

  float VB = D5*D2 - D1*D6;
  float w = D2/(D2-D6);
  float VA = D3*D6 - D5*D4;
  float x = (D4-D3)/((D4-D3)+(D5-D6));

  float VC = D1*D4 - D3*D2;
  float denom = 1.0f/(VA+VB+VC);
  point face1 = A + (VB*denom)*AB + (VC*denom)*AC;
  if (D1<= 0.0f && D2<= 0.0f)
    return A;
  
  else if(D3>=0.0f && D4<=D3)
    return B;
    
  else if(VC<=0.0f && D1>=0.0f && D3<=0.0f)
    return point(A+(v*AB));
      
  else if(D6>=0.0f && D5<=D6)
    return C;
        
  else if(VB<=0.0f && D2>=0.0f && D6<=0.0f)
    return point(A + w*AC);
              
  else if (VA<=0.0f && (D4-D3)>=0.0f && (D5-D6)>=0.0f)
    return point(B+ x*(C-B));
                  
  else
    return face1;
                      
                      
            

}
  
  
point closest_point_tetrahedron(const point &A, const point &B, const point &C, const point &D)
{
  vector Nabc = cross(B-A,C-A);
  vector Nacd = cross(C-A,D-A);
  vector Nadb = cross(D-A,B-A);
  vector AP = -A;
  vector AD = D-A;
  vector AB = B-A;
  vector AC = C-A;
  float min_dis = 20000;
  point clst_pnt{0,0,0};
                      

  if (dot(AP,Nabc)*dot(AD,Nabc)<=0.0f)
    {point temppoint = closest_point_triangle(A,B,C); /* point for face ABC*/
    float com_dist = norm(temppoint);
    if (com_dist<min_dis)
      {min_dis = com_dist;
      clst_pnt = temppoint;}}

  
  if (dot(AP,Nacd)*dot(AB,Nacd)<=0.0f)
    {point temppoint = closest_point_triangle(A,C,D); /* point for face ACD*/
    float com_dist = norm(temppoint);
    if (com_dist<min_dis)
      {min_dis = com_dist;
      clst_pnt = temppoint;}}
                        

  if (dot(AP,Nadb)*dot(AC,Nadb)<=0.0f)
    {point temppoint = closest_point_triangle(A,D,B); /* point for face ADB*/
    float com_dist = norm(temppoint);
    if (com_dist<min_dis)
      {min_dis = com_dist;
      clst_pnt = temppoint;}} 
  
  else {
    point temppoint = closest_point_triangle(B,D,C); /* point for face BDC*/
    float com_dist = norm(temppoint);
    if (com_dist<min_dis)
      {min_dis = com_dist;
      clst_pnt = temppoint;}}   

return clst_pnt;


}