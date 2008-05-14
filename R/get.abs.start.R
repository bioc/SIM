`get.abs.start` <-
function(chr, pos, dataset)
{
	if(is.numeric(chr))
	{
		msg = paste("the inserted chromosome column of the", dataset, "dataset is numeric")
		print(msg)
	}

	if(!is.numeric(chr))
	{
		print(paste("the inserted chromosome column of the", dataset, "is not numeric"))
		print("Trying to change x and y into numeric values")
		chr[chr == "X"] <- 24
		chr[chr == "Y"] <- 23
	
		if(sum(is.na(as.numeric(chr))) > 0)
		{	
			msg = paste("the inserted chromosome column of the", dataset, "dataset can not be converted, 
			please check your chromosome column and check for none integer values like NA's and remove them
			or convert them to an integer ")
			stop(msg)
		}
		chr <- as.numeric(chr)
		print("Convertion succesfull")
	}

	if(is.numeric(pos))
	{
		msg = paste("the inserted start position column of the", dataset, "dataset is numeric")
		print(msg)
	}

	if(!is.numeric(pos))
	{
		msg = paste("the inserted start postion column of the", dataset, "is not numeric,
		please check your position column and check for none integer values like NA's and remove them
		or convert them to an integer")
		stop(msg)
	}

	print("Trying to generate an absolute start object")


	Abs.start <- (chr * 1000000000) + pos
	print("Succesfully generated an absolute start object")
return(Abs.start)
}

