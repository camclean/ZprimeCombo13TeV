.text
.p2align 4,,
.globl atomic_add
.type atomic_add, @function
atomic_add:
    lock xaddq %rsi, (%rdi)
	ret
.size atomic_add, .-atomic_add

.p2align 4,,
.globl atomic_get
.type atomic_get, @function
atomic_get:
    movq (%rdi), %rax
	ret
.size atomic_get, .-atomic_get

.p2align 4,,
.globl atomic_set
.type atomic_set, @function
atomic_set:
    movq %rsi, (%rdi)
	ret
.size atomic_set, .-atomic_set

.p2align 4,,
.globl atomic_inc
.type atomic_inc, @function
atomic_inc:
    lock incq (%rdi)
    ret
.size atomic_inc, .-atomic_inc
