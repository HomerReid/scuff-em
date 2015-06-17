%s/<span class=codename>\([^<]*\)<\/span>/[[\1]]/g
%s/<span class=CodeName>\([^<]*\)<\/span>/[[\1]]/g
%s/<span class="CodeName">\([^<]*\)<\/span>/[[\1]]/g
%s/<code>\([^<]*\)<\/code>/``\1``/g
%s/<b>\([^<]*\)<\/b>/**\1**/g
%s/<i>\([^<]*\)<\/i>/*\1*/g
%s/<a class="CodeName" href="\([^"]*\)">\([^<]*\)<\/a>/[``\2''](\1)/g
%s/<a class="codename" href="\([^"]*\)">\([^<]*\)<\/a>/[``\2''](\1)/g
%s/^[ ]*//g
%s/^[ ]*<p>[ ]*$//g
%s/<pre class="Listing">/````/g
%s/<pre class=Listing>/````/g
%s/^<pre>/````/g
%s/<\/pre>/````/g
%s/&mu;/$\\mu$/g
%s/&eps;/$\\epsilon$/g
