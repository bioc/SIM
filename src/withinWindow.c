
void withinWindow(double *x, double *y, int *nx, int *ny, double *w)
{	
	int i, j, shift = 0;

	for(i = 0; i < *ny; i++)
	{
		j = shift;
		
		if(x[j] > y[i] + w[i + *ny]) //there are no x[j]'s left of y[i], move to next y[i] 
		{
			w[i] = w[i + *ny] = 0;
			continue;
		}
			
		while(x[j] < y[i] - w[i] && j < *nx)
		{
			if(i < *ny-1)           
			{
				if(x[j] < y[i+1] - w[i+1])
					shift = j; //keep track of a new start position for x[j]
			}
			j++;					
		}
	
		if(j == *nx)
		{
			w[i] = w[i + *ny] = 0;
		}
		else if(x[j] <=  y[i] + w[i + *ny]) 
		{
			w[i] = j + 1; //store left limit
			while(x[j] <= y[i] + w[i + *ny] && j < *nx)
				j++;				
			w[i + *ny] = j; //store right limit
		}
		else
			w[i] = w[i + *ny] = 0;
	}		
}



