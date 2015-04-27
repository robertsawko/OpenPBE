"Coding style 
"No tabs, just spaces (Google guide)
set expandtab number
"Indents are four columns (OF guide)
set shiftwidth=4
set softtabstop=4
"OF and Google guides
set textwidth=80

"N-s don't intent namespaces
set cino=>4,N-s,+4,(0,U1,W4,m1,l1,g1,h1,i4
set cinkeys=0{,0},0(,0),:,!^F,o,O,e

"Line length guard
highlight OverLength ctermbg=red ctermfg=white guibg=#592929
match OverLength /\%81v.\+/
set ruler number

"Indexed OF tags
set tags+=~/.vim/tags/of22


"wildignre allows us to ignore some files when searching and opening from
"within vim
set wildignore=*.o,*.obj,.git,lnInclude,linux64GccDPOpt,tags,*.dep,gmon.out,doc,build

"Compilation
set makeprg=cd\ build;make

set foldmethod=syntax

"A workaround for slow foldmethod
autocmd InsertEnter * let w:last_fdm=&foldmethod | setlocal foldmethod=manual
autocmd InsertLeave * let &l:foldmethod=w:last_fdm
