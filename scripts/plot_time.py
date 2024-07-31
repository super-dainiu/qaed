import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scienceplots


for name in ['skewrand']:
    table = pd.read_csv(f'{name}.csv')

    for loglog in [True, False]:
        with plt.style.context(['science', 'no-latex', 'ieee', 'grid', 'vibrant', 'high-vis', 'notebook']):
            fig, ax1 = plt.subplots(figsize=(8, 6))

            aed_table = table[table['deflation strategy'] == 'AED']
            iqr_table = table[table['deflation strategy'] == 'w/o AED']

            ax1.plot(aed_table['matrix size'], aed_table['total time'], label='AED Total Time', linewidth=3, marker='o', color='blue', linestyle='-')
            ax1.plot(iqr_table['matrix size'], iqr_table['total time'], label='IQR Total Time', linewidth=3, marker='o', color='green', linestyle='-')

            ax1.set_xlabel('Matrix Size', fontsize=20)
            ax1.set_ylabel('Time (s)', fontsize=20)
            if loglog:
                ax1.set_xscale('log')
                ax1.set_xticks([64, 128, 256, 512, 1024, 2048])
                ax1.set_xticklabels(['64', '128', '256', '512', '1024', '2048'], fontsize=12)
            else:
                ax1.set_xlim(0, 2100)
            ax1.set_yscale('log')
            ax1.legend(fontsize=18)

            if iqr_table['matrix size'].max() < aed_table['matrix size'].max():
                ax1.axvspan(iqr_table['matrix size'].max(), aed_table['matrix size'].max(), alpha=0.3, color='gray')

            fig.suptitle(f'{"Log-Log" if loglog else "Linear"} Plot of Execution Time', fontsize=24)

            plt.tight_layout()
            fname = f'{name}_time' + ('_loglog' if loglog else '') + '.pdf'
            fig.savefig(fname, bbox_inches='tight')
            plt.close(fig)

        with plt.style.context(['science', 'no-latex', 'ieee', 'grid', 'vibrant', 'high-vis', 'notebook']):
            fig, ax3 = plt.subplots(figsize=(8, 6))
            ax3.plot(aed_table['matrix size'], aed_table['total QR sweeps'], label='AED Total QR Sweeps', linewidth=3, linestyle='--', color='red', marker='*')
            ax3.plot(iqr_table['matrix size'], iqr_table['total QR sweeps'], label='IQR Total QR Sweeps', linewidth=3, linestyle='--', color='orange', marker='*')
            
            ax3.set_xlabel('Matrix Size', fontsize=20)
            ax3.set_ylabel('Total QR Sweeps', fontsize=20)
            if loglog:
                ax3.set_xscale('log')
                ax3.set_xticks([64, 128, 256, 512, 1024, 2048])
                ax3.set_xticklabels(['64', '128', '256', '512', '1024', '2048'], fontsize=12)
            else:
                ax3.set_xlim(0, 2100)
            ax3.set_yscale('linear')
            ax3.legend(fontsize=18)

            if iqr_table['matrix size'].max() < aed_table['matrix size'].max():
                ax3.axvspan(iqr_table['matrix size'].max(), aed_table['matrix size'].max(), alpha=0.3, color='gray')


            fig.suptitle(f'{"Log-Log" if loglog else "Linear"} Plot of QR Sweeps', fontsize=24)

            plt.tight_layout()
            fname = f'{name}_qr_sweeps' + ('_loglog' if loglog else '') + '.pdf'
            fig.savefig(fname, bbox_inches='tight')
            plt.close(fig)

        with plt.style.context(['science', 'no-latex', 'ieee', 'grid', 'vibrant', 'high-vis', 'notebook']):
            fig, ax1 = plt.subplots(figsize=(8, 6))

            # plot log-log of results (X: matrix size, Y: total time) vs (X: matrix size, Y: total time for w/o AED)
            aed_table = table[table['deflation strategy'] == 'AED']
            iqr_table = table[table['deflation strategy'] == 'w/o AED']

            ax1.plot(aed_table['matrix size'], aed_table['total time'], label='AED Total Time', linewidth=3, marker='o', color='blue', linestyle='-')
            ax1.plot(iqr_table['matrix size'], iqr_table['total time'], label='IQR Total Time', linewidth=3, marker='o', color='green', linestyle='-')

            ax1.plot(aed_table['matrix size'], aed_table['time to construct Q'], label='AED Time to Construct Q', linewidth=2, marker='s', color='blue', linestyle='-.')
            ax1.plot(iqr_table['matrix size'], iqr_table['time to construct Q'], label='IQR Time to Construct Q', linewidth=2, marker='s', color='green', linestyle='-.')

            ax1.plot(aed_table['matrix size'], aed_table['time for AED'], label='AED Time', linewidth=2, marker='^', color='blue', linestyle=':')

            if iqr_table['matrix size'].max() < aed_table['matrix size'].max():
                ax1.axvspan(iqr_table['matrix size'].max(), aed_table['matrix size'].max(), alpha=0.3, color='gray')

            ax1.set_xlabel('Matrix Size', fontsize=20)
            ax1.set_ylabel('Time (s)', fontsize=20)

            if loglog:
                ax1.set_xscale('log')
                ax1.set_xticks([64, 128, 256, 512, 1024, 2048])
                ax1.set_xticklabels(['64', '128', '256', '512', '1024', '2048'])
            else:
                ax1.set_xlim(0, 2100)

            ax1.set_yscale('log')
            ax1.set_ylim(1e-2, 1e6)
            ax1.tick_params(axis='y')

            ax2 = ax1.twinx()
            ax2.plot(aed_table['matrix size'], aed_table['total QR sweeps'], label='AED Total QR Sweeps', linewidth=3, linestyle='--', color='red', marker='*')
            ax2.plot(iqr_table['matrix size'], iqr_table['total QR sweeps'], label='IQR Total QR Sweeps', linewidth=3, linestyle='--', color='orange', marker='*')
            ax2.set_ylabel('Total QR Sweeps', fontsize=20)
            ax2.tick_params(axis='y')

            lines1, labels1 = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper left', bbox_to_anchor=(0.02, 0.98), fontsize=12)

            fig.suptitle('Serial Execution Time and Total QR Sweeps', fontsize=24)

            # save plot to file
            fname = f'{name}_time_sweeps' + ('_loglog' if loglog else '') + '.pdf'
            fig.savefig(fname, bbox_inches='tight')
            plt.close(fig)


    with plt.style.context(['science', 'no-latex', 'ieee', 'grid', 'vibrant', 'high-vis', 'notebook']):
        fig, ax1 = plt.subplots(figsize=(8, 6))

        aed_table = table[table['deflation strategy'] == 'AED']
        iqr_table = table[table['deflation strategy'] == 'w/o AED']

        matrix_sizes = np.sort(table['matrix size'].unique())
        tick_positions = range(len(matrix_sizes))

        normalized_aed_times = aed_table['total time'] / iqr_table.set_index('matrix size')['total time'].reindex(aed_table['matrix size']).values
        normalized_iqr_times = [1] * len(iqr_table)

        if iqr_table['matrix size'].max() < aed_table['matrix size'].max():
            normalized_iqr_times = np.append(normalized_iqr_times, [np.inf] * (len(aed_table) - len(iqr_table)))
            normalized_aed_times = np.append(normalized_aed_times[:len(iqr_table)], [normalized_aed_times[len(iqr_table)-1:len(iqr_table)]] * (len(aed_table) - len(iqr_table)))

        width = 0.35  # Bar width
        bar_aed = ax1.bar([pos - width/2 for pos in tick_positions], normalized_aed_times, width, label='AED Relative to IQR', color='blue')
        bar_iqr = ax1.bar([pos + width/2 for pos in tick_positions], normalized_iqr_times, width, label='IQR (baseline)', color='green')

        # add invalid label to iqr if it is np.inf
        for i, val in enumerate(bar_iqr):
            if val.get_height() == np.inf:
                # plot a grey bar with shadow
                ax1.bar(val.get_x() + val.get_width() / 2, 1, width, color='grey', alpha=0.3)
                ax1.text(val.get_x() + val.get_width() / 2, 0.2, 'Time Exceeded', ha='center', va='bottom', fontsize=12, color='green', rotation=90)

        ax1.set_xlabel('Matrix Size', fontsize=20)
        ax1.set_ylabel('Normalized Time', fontsize=20)
        fig.suptitle('Relative Execution Time Comparison', fontsize=24)
        ax1.legend(fontsize=18)

        ax1.set_xticks(tick_positions)
        ax1.set_xticklabels(matrix_sizes)  # Custom labels

        plt.tight_layout()
        fname = f'{name}_relative_time.pdf'
        fig.savefig(fname, bbox_inches='tight')
        plt.close(fig)