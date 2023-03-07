#include<iostream>
#include<cmath>

using namespace std;

int main()
{
   int enter_num, temp_num, sum = 0;
   int divisor, digit, count = 0;

   cout<<"Please enter number"<<"\n";
   cin>>enter_num;

   temp_num = enter_num;
//9 10
   // Counting the number of digits in the entered integer
   while (temp_num != 0)
   {
       temp_num = temp_num/10;
       cout << "temp_num: " << temp_num ;
       count++;
   }

   temp_num = enter_num;

   // Extracting the digits
   cout<<"Individual digits in the entered number are ";
   do
   {
       divisor = static_cast<int>(pow(10.0, --count));
       digit = temp_num / divisor;
       temp_num = temp_num % divisor;

       cout<<" "<<digit;
       sum = sum + digit;
   }
   while(count != 0);

   cout<<"\n"<<"Sum of the digits is = "<<sum<<"\n";

   return 0;
}