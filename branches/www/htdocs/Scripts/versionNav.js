function VersionNavigateur(Netscape, Explorer) 
  {
    if ((navigator.appVersion.substring(0,3) >= Netscape && navigator.appName == 'Netscape') ||      
        (navigator.appVersion.substring(0,3) >= Explorer && navigator.appName.substring(0,9) == 'Microsoft'))
           return true;
    else
	   return false;
  }
