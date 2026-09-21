.PHONY: figure validate clean

figure:
	python3 scripts/plot_figure2.py data/figure2_mach5_collision_production.csv -o output/figure2_mach5.png

validate:
	bash scripts/run_mach5.sh

clean:
	rm -f output/*.png output/*_regenerated.csv build/order8_m5_validation
	rm -rf runs
