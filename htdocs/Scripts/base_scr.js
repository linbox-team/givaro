if(document.images) {
        emailon = new Image();
        emailon.src = "Images/email2.gif"
        emailoff = new Image();
        emailoff.src = "Images/email.gif";
        emailline = new Image();
        emailline.src = "Images/emailline.gif";

        phoneon = new Image();
        phoneon.src = "Images/telephone2.gif"
        phoneoff = new Image();
        phoneoff.src = "Images/telephonex.gif";
        phoneline = new Image();
        phoneline.src = "Images/phoneline.gif";

        faxon = new Image();
        faxon.src = "Images/fax2.gif"
        faxoff = new Image();
        faxoff.src = "Images/fax.gif";
        faxline = new Image();
        faxline.src = "Images/faxline.gif";

        adresson = new Image();
        adresson.src = "Images/adress2.gif"
        adressoff = new Image();
        adressoff.src = "Images/adress.gif";
        adressline = new Image();
        adressline.src = "Images/adressline.gif";
}


function msoverEMAIL(imgName) {
        if (document.images) {
                imgOn = eval(imgName + "on.src");
                document [imgName].src = imgOn;
        }
        cipher = "Ojfs-Lznqqfzrj.Izrfx@nrfl.kw"
        alpha = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
        key = 5
        clear = ""
        for (i = 0; i < cipher.length; i++) {
            if (alpha.indexOf (cipher.charAt (i)) == -1) {
                clear += cipher.charAt (i)
            } else {     
                offset = (alpha.indexOf (cipher.charAt (i)) + alpha.length - key) % alpha.length
                clear += alpha.charAt (offset)
            }                               
        }
        affiche('<font face="Arial,Helvetica"><font size="+0"><big>' + clear + '</big></font></font>');
}

function msoverPHONE(imgName) {
        if (document.images) {
                imgOn = eval(imgName + "on.src");
                document [imgName].src = imgOn;
        }
        affiche('<font face="Arial,Helvetica"><font size="+0"><big>+33 (0) 476 514 866</big></font></font>');
}

function msoverFAX(imgName) {
        if (document.images) {
                imgOn = eval(imgName + "on.src");
                document [imgName].src = imgOn;
        }
        affiche('<font face="Arial,Helvetica"><font size="+0"><big>+33 (0) 476 631 263</big></font></font>');
}

function msoverADD(imgName) {
        if (document.images) {
                imgOn = eval(imgName + "on.src");
                document [imgName].src = imgOn;
        }
        affiche('<font face="Arial,Helvetica"><font size="+0"><big>Laboratoire Jean Kuntzmann<br>51, av. des Mathématiques<br>BP 53X, 38041 Grenoble Cedex 9<br>FRANCE</big></font></font>');
}

function msover(imgName,text) {
        if (document.images) {
                imgOn = eval(imgName + "on.src");
                document [imgName].src = imgOn;
                lineOn = eval(imgName + "line.src");
                document.line.src = lineOn;
        }
        affiche(text);
}

function msout(imgName) {
        if (document.images) {
                imgOff = eval(imgName + "off.src");
                document [imgName].src = imgOff;
        }
        cache(' ');
}      

function affiche(text) {
    var object = document.getElementById('Line');
    object.innerHTML=text;
}
      
function cache() {
    var object = document.getElementById('Line');
    object.innerHTML=" ";
}


// Anti-Spam Script Adapted from Laurent Fousse (c) 2007
function CryptMail(AddOn) {
    cipher = "<f mwjk=\"rfnqyt:Ojfs-Lznqqfzrj.Izrfx@nrfl.kw\" " + AddOn + ">Ojfs-Lznqqfzrj.Izrfx@nrfl.kw</f>"
    alpha = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    key = 5
    clear = ""
    for (i = 0; i < cipher.length; i++) {
	if (alpha.indexOf (cipher.charAt (i)) == -1) {
	    clear += cipher.charAt (i)
	} else {     
	    offset = (alpha.indexOf (cipher.charAt (i)) + alpha.length - key) % alpha.length
	    clear += alpha.charAt (offset)
	}                               
    }
    document.write (clear)
}

function emailtomecrypted() {
    cipher = "rfnqyt:Ojfs-Lznqqfzrj.Izrfx@nrfl.kw"
    alpha = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    key = 5
    clear = ""
    for (i = 0; i < cipher.length; i++) {
	if (alpha.indexOf (cipher.charAt (i)) == -1) {
	    clear += cipher.charAt (i)
	} else {     
	    offset = (alpha.indexOf (cipher.charAt (i)) + alpha.length - key) % alpha.length
	    clear += alpha.charAt (offset)
	}                               
    }
    window.location.href = clear;
}
