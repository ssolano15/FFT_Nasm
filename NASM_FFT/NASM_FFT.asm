extern calloc ; Importar la funcion Calloc
extern free ; Importar la funcion free

calloc_int_size: ;def funcion
	push 4 ;sizeof(int)
	push qword [rbp - 56] ;tamaño
	call calloc ;llama a calloc
	add rsp, 16 ;pop que no guarda datos
	ret ;return

calloc_double_2size: ;def funcion
	push 8 ;sizeof(double) 
	push qword [rbp - 56] ;tamaño
	shl qword [rsp], 1 ;multiplica tamaño por dos
	call calloc ;llama a calloc
	add rsp, 16 ;pop que no guarda datos
	ret ;return

section   .data ; No Progra
				message:  db        "Terminó FFT", 10 ; No Progra


section .text
	
	global fft ;define ffft

fft:
	;Guardar parametros in_data y size
	push   rbp
	mov    rbp,rsp
	sub    rsp, 48
	push   rbx
	push   rsi
	push   rdi

	xor rcx, rcx; int k=0
	dec rcx
	log_loop: ; while ((1 << k) < size)_
		inc rcx ;suma al log_loop

		; (1 << k) < size
		xor rax, rax
		inc rax
		shl rax, cl

		cmp rax, [rbp - 56] ;compara eax con el tamaño
		jb log_loop	; si es menor salta

	mov qword[rbp-8], rcx; se asignan los primeros -8 para k

	; int *rev = (int*) calloc(size, sizeof(int))
	call calloc_int_size
	mov [rbp - 16], rax

	mov rdx, [rbp - 16]
	mov qword [rdx], 0; rev[0]=0
	mov qword [rbp - 24], -1; int high1=-1 *************************

	xor rcx, rcx
	inc rcx
	rev_loop: ;for (i = 1; i < size; ++i)		
		; if (i & (i - 1) == 0)
		mov rax, rcx
		dec rax
		test rax, rcx
		jnz .false
			inc qword [rbp - 24]
		.false:

		;i ^ (1 << high1)
		xor rax, rax
		inc rax 
		push rcx
		mov rcx, qword [rbp - 24]
		shl rax, cl
		pop rcx
		xor rax, rcx

		mov rdx, [rbp - 16]
		mov rbx, [rdx + 8 * rax]; ****************
		mov qword [rdx + 8 * rcx], rbx ;*******************

		;1<<(k-high1-1)
		mov rdx, qword [rbp - 8]
		sub rdx, qword [rbp - 24]
		dec rdx
		xor rax, rax
		inc rax
		push rcx
		mov rcx, rdx
		shl rax, cl
		pop rcx

		mov rdx, [rbp - 16]
		xor qword[rdx + 8 * rcx], rax;*********************

		inc rcx
		cmp rcx, [rbp - 56]
		jb rev_loop

	;double *roots = (double*) calloc(2 * size, sizeof(double))
	call calloc_double_2size
	mov [rbp-32], rax

	;double alpha= 2*M_PI/ size
	fldpi
	push 2
	fmul qword[rsp]
	add rsp, 8
	fdiv qword[rbp - 56]

	; for (i = 0; i < size; ++i)
	xor rcx, rcx
	roots_loop:
		fld st0
		push rcx
		fmul qword [rsp]
		pop rcx
		fsincos

		mov rax, rcx
		shl rax, 1
		mov rdx, [rbp - 32]
		fstp qword [rdx + 8 * rax];***********
		fstp qword [rdx + 8 * rax + 8];************

		inc rcx
		cmp rcx, [rbp - 56]
		jb roots_loop

	; double *cur = (double*) calloc(2 * size, sizeof(double));
	call calloc_double_2size
	mov [rbp-40], rax

	;for (i = 0; i < size; ++i)
	xor rcx, rcx
	cur_loop:
		;int ni = rev[i]
		mov rdx, [rbp - 16]
		mov rbx, [rdx + 4 * rcx]

		shl rbx, 1
		mov rdx, [rbp - 64]
		movsd xmm0, [rdx + 8 * rbx];***********
		movsd xmm1, [rdx + 8 * rbx + 8]

		shl rcx, 1
		mov rdx, [rbp - 40]
		movsd [rdx + 8 * rcx], xmm0
		movsd [rdx + 8 * rcx], xmm1
		shr rcx, 1

		inc rcx
		cmp rcx, [rbp - 56]
		jb cur_loop

	;free (rev)
	push qword [rbp-16]
	call free
	add rsp, 8

	; for (len = 1; len < size; len <<= 1)
	; for (len = 1; len < size; len <<= 1)
	xor rcx, rcx
	inc rcx
	len_loop:
		; double *ncur = (double*) calloc(2 * size, sizeof(double))
		push rcx
		call calloc_double_2size
		mov [rbp - 48], rax
		pop rcx
		
		mov rax, qword [rbp -56]
		div rcx
		shr rax, 1
		
		xor rdx, rdx
		p1_loop:
			xor rbx, rbx
			i_loop:
				mov rdi, [rbp - 32]
				mov rsi, rbx
				imul rsi, rax
				shl rsi, 1
				
				movsd xmm2, [rdi + 8 * rsi];***********
				movsd xmm3, [rdi + 8 * rsi + 8];************
				
				mov rdi, [rbp - 40]
				mov rsi, rdx
				add rsi, rcx
				shl rsi, 1
				
				movsd xmm4, [rdi + 8 * rsi]
				movsd xmm5, [rdi + 8 * rsi + 8]
				
				movsd xmm6, xmm2
				mulsd xmm6, xmm4
				movsd xmm7, xmm3
				mulsd xmm7, xmm5
				
				movsd xmm0, xmm6
				subsd xmm0, xmm7
				
				movsd xmm6, xmm2
				mulsd xmm6, xmm5
				movsd xmm7, xmm3
				mulsd xmm7, xmm4
				
				movsd xmm1, xmm6
				addsd xmm1, xmm7
				
				mov rdi, [rbp - 40]
				mov rsi, rdx
				shl rsi, 1
				
				movsd xmm2, [rdi + 8 * rsi]
				movsd xmm3, [rdi + 8 * rsi + 8]
				
				mov rdi, [rbp - 48]
				movsd xmm4, xmm2
				addsd xmm4, xmm0
				movsd xmm5, xmm3
				addsd xmm5, xmm1
				
				movlpd [rdi + 8 * rsi], xmm4
				movlpd [rdi + 8 * rsi + 8], xmm5
				
				shr rsi, 1
				add rsi, rcx
				shl rsi, 1
				
				movsd xmm6, xmm2
				subsd xmm6, xmm0
				movsd xmm7, xmm3
				subsd xmm7, xmm1
				
				movlpd [rdi + 8 * rsi], xmm6
				movlpd [rdi + 8 * rsi + 8], xmm7
								
				inc rbx
				inc rdx
				cmp rbx, rcx
				jb i_loop
			
			
			add rdx, rcx
			cmp rdx, [rbp - 56]
			jb p1_loop
		
		; free(ncur)
		push rcx
		push qword [rbp - 40]
		call free
		add rsp, 8
		pop rcx
		
		mov rdx, [rbp - 48]
		mov [rbp - 40], rdx
		
		shl rcx, 1
		cmp rcx, [rbp - 56]
		jb len_loop

	; free(roots)
	push qword [rbp - 32]
	call free
	add rsp, 8

	; return cur;
	mov rax, [rbp - 40]

	pop rdi
	pop rsi
	pop rbx
	mov rsp, rbp
	pop rbp
	ret	
   	syscall                           ; invoke operating system to do the write      No Progra
    	mov       rax, 60                 ; system call for exit      No Progra
    	xor       rdi, rdi                ; exit code 0      No Progra
    	syscall 
    
    
    
    ; rbp -> almacenar datos
    ; ecx -> controlar loops
