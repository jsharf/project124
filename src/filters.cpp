#include <vector>
using namespace std;

// Update each position as the average of its neighbor values
void blurFilter(vector <int>& list)
{
    int last = 0;
    for (int i=0; i<list.size(); i++)
    {
        int tmp = list[i];
        list[i] = (list[i+1] + list[i] + last)/3;
        last = tmp;
    }
}

// takes discrete derivative of each position in the input list
// for each value list[i]: list[i] = list[i] - list[i-1]
void sobel(vector <int>& list)
{
    int last = 0;
    for (int i=0; i<list.size(); i++)
    {
        int tmp = list[i];
        list[i] = list[i] - last;
        last = tmp;
    }
}
