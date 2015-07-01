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

