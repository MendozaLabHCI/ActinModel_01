function FormatedNumber = SimFormat(Num)
    
    FormatedNumber = sprintf('%#0.3f',Num);
    FormatedNumber(FormatedNumber == '.') = 'd';
    

end