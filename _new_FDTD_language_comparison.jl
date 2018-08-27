using Printf

const NX = 300::Int64								# 空間セル数 X [pixels]
const NY = 400::Int64								# 空間セル数 Y [pixels]

const dx = 0.01::Float64								# 空間刻み [m]
const dt = 20.0e-6::Float64							# 時間刻み [s]

const Nstep = 10000::Int64								# 計算ステップ数 [回]

const freq = 1.0e3::Float64							# 初期波形の周波数 [Hz]

const ρ = 1.3::Float64								# 密度ρ [kg/m^3]
const κ = 142.0e3::Float64							# 体積弾性率κ [Pa]



# 事前準備 #########################################################
waveformfile = open("waveform.txt", "w")

function calc()
# メインループ #########################################################
    Vx = zeros(Float64, NX+1,NY  )			# x方向粒子速度 [m/s]
    Vy = zeros(Float64, NX,  NY+1)			# y方向粒子速度 [m/s]
    P  = zeros(Float64, NX,  NY  )			# 音圧 [Pa]
for n = 0:Nstep

	# 更新（ここが FDTD の本体）
	# 粒子速度の更新
    @inbounds for i=2:NX, j=1:NY
        Vx[i,j] += - dt / (ρ * dx) * ( P[i,j] - P[i-1,j] );
    end
    @inbounds for i=1:NX, j=2:NY
        Vy[i,j] += - dt / (ρ * dx) * ( P[i,j] - P[i,j-1] );
	end
    # 音圧の更新
    @inbounds for i=1:NX, j=1:NY
	    P[i,j] += - ( κ * dt / dx ) * ( ( Vx[i+1,j] - Vx[i,j] ) + ( Vy[i,j+1] - Vy[i,j] ) );
    end
	# 初期波形を準備（正弦波×１波 with ハン窓）
	if n < (1.0/freq)/dt
		sig = (1.0-cos(2.0*pi*freq*n*dt))/2.0 * sin(2.0*pi*freq*n*dt)
	else
		sig = 0.0
	end

	# 音源
	@inbounds P[fld(NX,4)+1,fld(NY,3)+1] = sig
	@printf("%5d / %5d\r", n, Nstep);
#=
	# 波形ファイル出力（時刻, 音源, 中央点の音圧）
	write(waveformfile,"$(dt*n)\t$sig\t$(P[Int32(floor(NX/2+1)),Int32(floor(NY/2+1))])\n")

	# 音圧分布ファイル出力（50ステップ毎）
	if n % 50 == 0
		fieldfilename = @sprintf("field%06d.txt",n)
		fieldfile = open(fieldfilename,"w")
		for i=1:NX
			for j=1:NY
				write(fieldfile,"$(P[i,j])\t")
			end
			write(fieldfile,"\n")
		end
		close(fieldfile)
	end
=#
end
end

calc()
# 事後処理 #########################################################
close(waveformfile)
