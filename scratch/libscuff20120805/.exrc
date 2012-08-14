if &cp | set nocp | endif
map  :cn 
map  "_ddP
map  :!remake
map  :tab split:exec("tag ".expand("<cword>"))
let s:cpo_save=&cpo
set cpo&vim
nmap gx <Plug>NetrwBrowseX
map nf :e <cfile> 
map <F1> :w  :make  
nnoremap <silent> <Plug>NetrwBrowseX :call netrw#NetrwBrowseX(expand("<cWORD>"),0)
map <F5> :cn 
map <F4> :w  :!./m                
map <F3> :.,/end.align/s/^%//g 
map <F2> :.,/end.align/s/^/%/g 
map Ý :vsp :exec("tag ".expand("<cword>"))
map Û :sp :exec("tag ".expand("<cword>"))
map ì :CN
map ë 50k
map ê 50j
cmap dddf .,/^}/d
cmap yyyf .,/^}/y
let &cpo=s:cpo_save
unlet s:cpo_save
set background=dark
set fileencodings=ucs-bom,utf-8,default,latin1
set guifont=Courier\ 10\ Pitch\ 12
set guioptions=aegimrLt
set helplang=en
set iconstring=vim:\ %t
set mouse=a
set shellcmdflag=-lc
set splitright
set termencoding=utf-8
set visualbell
set window=56
" vim: set ft=vim :
