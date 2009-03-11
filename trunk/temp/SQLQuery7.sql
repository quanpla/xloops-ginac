USE [testDB]
GO
/****** Object:  UserDefinedFunction [dbo].[fn_ConvertInt2Char_RND]    Script Date: 03/11/2009 17:17:41 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE FUNCTION [dbo].[fn_ConvertInt2Char_RND](@IntNum INT)
RETURNS CHAR(1)
AS BEGIN
/*************************************************************
**	int		0	->	9	10	->	35	36	->	61
**	char	0	...	9	A	...	Z	a	...	z
**	ascii	48	...	57	65	...	90	97	...	122
*************************************************************/
	DECLARE @ConvertInt INT;
	IF @IntNum > 61
	BEGIN
		SET @ConvertInt = NULL;
	END
	ELSE
		IF @IntNum < 10
		BEGIN
			SET @ConvertInt = @IntNum + 48;
		END -- is a number
		ELSE
		BEGIN -- not a number
			IF @intNum < 36
			BEGIN	--	is CAP
				SET @ConvertInt = @IntNum + 65 - 10;
			END
			ELSE
			BEGIN	--	not CAP
				SET @ConvertInt = @IntNum + 97 - 36;
			END
		END

	RETURN CHAR(@ConvertInt);
END
GO

USE [testDB]
GO
/****** Object:  UserDefinedFunction [dbo].[rand_Str]    Script Date: 03/11/2009 17:17:54 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE FUNCTION [dbo].[rand_Str](@strMaxLen INT, @isFixLen INT = 0)
RETURNS VARCHAR(2000)
AS BEGIN
/***********************************************************
** Generate a string with max len input
**
**	SELECT dbo.rand_Str(12)
***********************************************
**	Input:
**	@strMaxLen
**	@isFixLen: if = 0: the return will have random len up to max len
				  = 1: the return string is strictly Max Len
***********************************************************/

--	The Return string
	DECLARE @return_data VARCHAR(2000), @return_length INT;
	SET @return_data = '';

	SELECT @return_length = CASE @isFixLen	WHEN 1 THEN @strMaxLen
											WHEN 0 THEN 1 + RandomNumber * (@strMaxLen - 1)
							END
	FROM dbo.rand_Number;

--	The main loop to generate data
	DECLARE @randInt INT;

	WHILE LEN(RTRIM(@return_data)) < @return_length
	BEGIN
		-- Generate a random ASCII, remember that rand function can not be called within one function
		--	We have to use a trick to use the view
		SELECT @randInt = FLOOR(RandomNumber*61.5)
		FROM dbo.rand_Number; -- this is the random number view

		SET @return_data = @return_data + ISNULL(dbo.fn_ConvertInt2Char_RND(@randInt), '');

	END	-- while < max len

	RETURN @return_data;
END
