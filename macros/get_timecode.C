/*This macro is used to convert pedVStime x-coordinate id into*/ 
/*human-readable date-time.                                   */
/*It also has two utility functions to help timecode-date converting*/

UInt_t convert_dateTotimecode(const char* date)//example:20140410162029
{
    if(strlen(date) != 14){
        printf("error: invalid date format!\n");
        printf("The expected format is: 20140410162029\n");
        return -1;
    }

    UInt_t timecode;
    int leapyear_num;
    int totalyear_num;
    char buffer[10];
    int year,init_year;
    int month,init_month;
    int day,init_day;
    int hour,init_hour;
    int miniute,init_miniute;
    int second,init_second;

    init_year=2013;
    init_month=1;
    init_day=1;
    init_hour=0;
    init_miniute=0;
    init_second=0;
    int monthday[12]={31,28,31,30,31,30,31,31,30,31,30,31};
    int leap_monthday[12]={31,29,31,30,31,30,31,31,30,31,30,31};
    int leap_flag=0;
    //-------
    strncpy(buffer,date,4);
    buffer[4]='\0';
    year=atoi(buffer);

    strncpy(buffer,date+4,2);
    buffer[2]='\0';
    month=atoi(buffer);

    strncpy(buffer,date+6,2);
    buffer[2]='\0';
    day=atoi(buffer);

    strncpy(buffer,date+8,2);
    buffer[2]='\0';
    hour=atoi(buffer);

    strncpy(buffer,date+10,2);
    buffer[2]='\0';
    miniute=atoi(buffer);

    strncpy(buffer,date+12,2);
    buffer[2]='\0';
    second=atoi(buffer);
    //printf("input:%s\n",date);
    //printf("output:%4.d%2.d%2.d%2.d%2.d%2.d\n",year,month,day,hour,miniute,second);
    //--------
    totalyear_num=year-init_year;
    leapyear_num=0;
    int current_year;
    for(int i=0;i<totalyear_num;i++){
        current_year=init_year+i;
        if(((current_year%4)== 0)&&((current_year%100)>0)){
            leapyear_num++;
        }
        else if((current_year%400)==0){
            leapyear_num++;
        }
    }
    if(((year%4)== 0)&&((year%100)!=0)){
        leap_flag=1;
    }
    else if((year%400)==0){
        leap_flag=1;
    }
    //----------
    int month_num=month-init_month;
    int day_num=day-1;
    timecode=leapyear_num*366*24*60*60 + (totalyear_num-leapyear_num)*365*24*60*60;
    if(leap_flag){
        for(int i=0;i<month_num;i++){
            day_num+=leap_monthday[i];
        }
    }
    else{
        for(int i=0;i<month_num;i++){
            day_num+=monthday[i];
        }
    }
    timecode+=day_num*24*60*60 + hour*60*60 + miniute*60 + second;

    return timecode;
}

bool check_leapyear(int year)
{
    if((year%4==0)&&(year%100)!=0){
        return true;
    }
    else if((year%400)==0){
        return true;
    }
    else{
        return false;
    }
}

char* convert_timecodeTodate(UInt_t timecode)
{
    static char date[50];

    int totalyear_num;
    int current_year;
    int current_month;
    int current_day;
    int current_hour;
    int current_miniute;
    int current_second;
    bool currentyear_leapflag=false;
    UInt_t totalsecond;
    int init_year=2013;
    int normal_monthday[12]={31,28,31,30,31,30,31,31,30,31,30,31};
    int leap_monthday[12]={31,29,31,30,31,30,31,31,30,31,30,31};
    int noramlyear_second=365*24*3600;
    int leapyear_second=366*24*3600;

    int temp_year;
    totalyear_num=-1;
    totalsecond=0;
    int tmp_totalsecond;
    while(timecode>=totalsecond){
        totalyear_num++;
        temp_year=init_year+totalyear_num;
        tmp_totalsecond=totalsecond;
        if(check_leapyear(temp_year)){
            totalsecond+=leapyear_second;
        }
        else{
            totalsecond+=noramlyear_second;
        }
    }
    current_year=temp_year;

    if(check_leapyear(current_year)){
        currentyear_leapflag=true;
    }
    int current_daynum=(timecode-tmp_totalsecond)/(24*3600);
    current_month=0;
    int temp_daynum=0;
    int temp_daynum_before;
    for(int i=0;i<12;i++){
        current_month++;
        if(currentyear_leapflag){
            temp_daynum+=leap_monthday[i];
        }
        else{
            temp_daynum+=normal_monthday[i];
        }
        if(temp_daynum > current_daynum){
            current_day=current_daynum-temp_daynum_before+1;
            break;
        }
        temp_daynum_before=temp_daynum;
    }

    int remaining_second=(timecode-tmp_totalsecond)%(24*3600);
    current_hour=remaining_second/3600;
    remaining_second-=current_hour*3600;
    current_miniute=remaining_second/60;
    current_second=remaining_second%60;

    sprintf(date,"%d/%d/%d/%d/%d/%d",current_year,current_month,current_day,current_hour,current_miniute,current_second);

    return date;
}


int get_timecode(char* filename,int seg_id)
{
    TFile f(filename);
    UInt_t time_second;
    TTree* tree=(TTree*)f.Get("PSD");
    tree->SetBranchAddress("time_second",&time_second);

    tree->GetEntry(seg_id*2000);
    printf("seg_id=%d, corresponding to %u seconds from the initial time\n",seg_id,time_second);
	printf("Date is (from 2013/01/01,00:00:00): %s\n",convert_timecodeTodate(time_second));
	
	return 0;
}
