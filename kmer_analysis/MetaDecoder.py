import subprocess


class MetaDecoder:
    def __init__(self, thread_count=4):

        self.thread_count = thread_count
        # self.fasta_path = fasta_path

    def print_subprocess(self, completed_SP: subprocess.CompletedProcess) -> bool:
        print(f'Subprocess\n {"=" * 10}\n')
        print(f'Exit Code (init): {completed_SP.returncode}')
        print(f'STDOUT (init): {completed_SP.stdout}')

        err = completed_SP.stderr
        if err is not None:

            print(f'STDERR (init): {completed_SP.stderr}')
            return True
        else:
            return False

    def write_data(self):
        pass

    def get_data(self):
        """
        Gets data from metadecoder file and puts them in data structures
        :return:
        """

        pass

    def cluster(self, fasta_path):
        output_prefix = fasta_path.split('.')[0]
        output_seed = "".join([output_prefix, '.SEED'])
        output_cluster = "".join([output_prefix, '.COVERAGE'])

        seed_cmd = ['metaencoder', 'seed', '--threads', self.thread_count,
                    '-f', fasta_path, '-o', ]

        seed_process = subprocess.run(seed_cmd)
        self.print_subprocess(completed_SP=seed_process)

        cluster_cmd = ['metadecoder', 'cluster', '-f', fasta_path, '-c', output_cluster, '-s', output_seed, '-o',
                       'METADECODER']

        cluster_process = subprocess.run(cluster_cmd)

        self.print_subprocess(completed_SP=cluster_process)

        # metadecoder seed --threads 50 -f ASSEMBLY.FASTA -o METADECODER.SEED
        # metadecoder cluster -f ASSEMBLY.FASTA -c METADECODER.COVERAGE -s METADECODER.SEED -o METADECODER
        return seed_process, cluster_process

if __name__ == '__main__':
    new_path = 'assembly.fasta'
    MD = MetaDecoder(thread_count=8)
    (seed, run) = MD.cluster(new_path)

    MD.print_subprocess(MD.completed_process)
