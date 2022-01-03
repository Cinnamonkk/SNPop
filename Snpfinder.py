import sys
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QDir, pyqtSlot, pyqtSignal
import requests
import pandas as pd
import json
import traceback


# Worker class signal handler
###########################################################################################
class WorkerSignals(QtCore.QObject):
    '''Signal Handler for the Worker Thread'''

    started = pyqtSignal()
    finished = pyqtSignal()
    result = pyqtSignal(object)
    error = pyqtSignal(tuple)

# Workers class/Thread
############################################################################################


class Worker(QtCore.QRunnable):
    '''Allows the tasks passed in to run in a seperate thread to prevent GUI freezing'''

    def __init__(self, function, *args, **kwargs):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)

        self.function = function
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''

        # Retrieve args/kwargs here and fire processing using them

        self.signals.started.emit()
        try:
            result = self.function(
                *self.args, **self.kwargs
            )
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            # Return the result of the processing
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit()  # Done


class Ui_MainWindow(QtWidgets.QMainWindow):

    # Worker thread for the UID summary request
    #################################################

    def UIDclicked(self, item):
        '''connects the thread to the infosum function'''
        worker = Worker(self.infosum, item)
        worker.signals.started.connect(self.infosearch)
        self.threadpool.start(worker)

        # Worker thread for the UID list request
        ##############################################################################################################################

    def RetrieveSNPs_clicked(self):
        '''When the Retrieve UIDs button is clicked empty the SNP list and load a new one using the available_SNV function'''
        self.UIDlist.clear()
        if self.Retmax.value() == 0:
            msg = QtWidgets.QMessageBox.information(self, ' ', 'Retmax value cant be 0,\nplease enter the number of results you want.',
                                                    QtWidgets.QMessageBox.Yes)

        else:
            worker = Worker(self.available_SNV,
                            self.Retstart.value(), self.Retmax.value(), self.genename.text(), self.clinsignificance.currentText())
            worker.signals.started.connect(self.Ui_Load_UIDs_off)
            worker.signals.result.connect(self.UIdlistfunc)
            worker.signals.finished.connect(self.Ui_Load_UIDs_on)
            self.threadpool.start(worker)

    def UIdlistfunc(self, Uidlist):
        '''Updates the UID list widget'''
        strlist = [str(x) for x in Uidlist]
        self.UIDlist.clear()
        self.UIDlist.addItems(tuple(strlist))


# Worker thread for the retrieve the frequencies request
#######################################################################


    def Get_Data_Clicked(self):
        '''Retrieves data for the selected rsIDs'''
        if self.SelectedUID.count() == 0:
            msg = QtWidgets.QMessageBox.warning(self, 'ERROR', 'Zero UIDs selected',
                                                QtWidgets.QMessageBox.Ok)
        else:
            self.DataTable.model().removeRows(0, self.DataTable.rowCount())
            worker1 = Worker(self.Popfinder, self.Selected_UID_list())
            worker2 = Worker(self.Phenfinder, self.Selected_UID_list())
            worker3 = Worker(self.Genotypefinder, self.Selected_UID_list())
            worker2.signals.started.connect(self.Ui_Get_Freqs_off)
            worker3.signals.started.connect(self.Ui_Get_Freqs_off)
            worker1.signals.started.connect(self.Ui_Get_Freqs_off)
            worker1.signals.result.connect(self.Tablemaker)
            worker2.signals.result.connect(self.Tablemaker2)
            worker3.signals.result.connect(self.Tablemaker3)
            worker1.signals.finished.connect(self.Ui_Get_Freqs_on)
            self.threadpool.start(worker1)
            self.threadpool.start(worker2)
            self.threadpool.start(worker3)

    def Tablemaker(self, alistoflists):
        '''Passes the SNP data into the Table'''
        self.DataTable.setRowCount(len(alistoflists[0]))
        for row in range(len(alistoflists[0])):
            self.DataTable.setItem(
                row, 0, QtWidgets.QTableWidgetItem(alistoflists[0][row]))
            self.DataTable.setItem(
                row, 1, QtWidgets.QTableWidgetItem(alistoflists[1][row]))
            self.DataTable.setItem(
                row, 2, QtWidgets.QTableWidgetItem(alistoflists[2][row]))
            self.DataTable.setItem(
                row, 3, QtWidgets.QTableWidgetItem(alistoflists[3][row]))
            self.DataTable.setItem(
                row, 4, QtWidgets.QTableWidgetItem(alistoflists[4][row]))
            try:
                self.DataTable.setItem(
                    row, 5, QtWidgets.QTableWidgetItem(str(alistoflists[5][row]['1000GENOMES:phase_3:ALL_minor'])))
            except KeyError:
                self.DataTable.setItem(
                    row, 5, QtWidgets.QTableWidgetItem('NA'))
            try:
                self.DataTable.setItem(
                    row, 6, QtWidgets.QTableWidgetItem(str(alistoflists[5][row]['1000GENOMES:phase_3:ALL_major'])))
            except KeyError:
                self.DataTable.setItem(
                    row, 6, QtWidgets.QTableWidgetItem('NA'))
            try:
                self.DataTable.setItem(
                    row, 7, QtWidgets.QTableWidgetItem(str(alistoflists[5][row]['1000GENOMES:phase_3:AFR_minor'])))
            except KeyError:
                self.DataTable.setItem(
                    row, 7, QtWidgets.QTableWidgetItem('NA'))
            try:
                self.DataTable.setItem(
                    row, 8, QtWidgets.QTableWidgetItem(str(alistoflists[5][row]['1000GENOMES:phase_3:AFR_major'])))
            except KeyError:
                self.DataTable.setItem(
                    row, 8, QtWidgets.QTableWidgetItem('NA'))
            try:
                self.DataTable.setItem(
                    row, 9, QtWidgets.QTableWidgetItem(str(alistoflists[5][row]['1000GENOMES:phase_3:EUR_minor'])))
            except KeyError:
                self.DataTable.setItem(
                    row, 9, QtWidgets.QTableWidgetItem('NA'))
            try:
                self.DataTable.setItem(
                    row, 10, QtWidgets.QTableWidgetItem(str(alistoflists[5][row]['1000GENOMES:phase_3:EUR_major'])))
            except KeyError:
                self.DataTable.setItem(
                    row, 10, QtWidgets.QTableWidgetItem('NA'))
            try:
                self.DataTable.setItem(
                    row, 11, QtWidgets.QTableWidgetItem(str(alistoflists[5][row]['1000GENOMES:phase_3:AMR_minor'])))
            except KeyError:
                self.DataTable.setItem(
                    row, 11, QtWidgets.QTableWidgetItem('NA'))
            try:
                self.DataTable.setItem(
                    row, 12, QtWidgets.QTableWidgetItem(str(alistoflists[5][row]['1000GENOMES:phase_3:AMR_major'])))
            except KeyError:
                self.DataTable.setItem(
                    row, 12, QtWidgets.QTableWidgetItem('NA'))
            try:
                self.DataTable.setItem(
                    row, 13, QtWidgets.QTableWidgetItem(str(alistoflists[5][row]['1000GENOMES:phase_3:EAS_minor'])))
            except KeyError:
                self.DataTable.setItem(
                    row, 13, QtWidgets.QTableWidgetItem('NA'))
            try:
                self.DataTable.setItem(
                    row, 14, QtWidgets.QTableWidgetItem(str(alistoflists[5][row]['1000GENOMES:phase_3:EAS_major'])))
            except KeyError:
                self.DataTable.setItem(
                    row, 14, QtWidgets.QTableWidgetItem('NA'))
            try:
                self.DataTable.setItem(
                    row, 15, QtWidgets.QTableWidgetItem(str(alistoflists[5][row]['1000GENOMES:phase_3:SAS_minor'])))
            except KeyError:
                self.DataTable.setItem(
                    row, 15, QtWidgets.QTableWidgetItem('NA'))
            try:
                self.DataTable.setItem(
                    row, 16, QtWidgets.QTableWidgetItem(str(alistoflists[5][row]['1000GENOMES:phase_3:SAS_major'])))
            except KeyError:
                self.DataTable.setItem(
                    row, 16, QtWidgets.QTableWidgetItem('NA'))

    def Tablemaker2(self, alistoflists):
        self.DataTable.setRowCount(len(alistoflists[0]))
        for row in range(len(alistoflists[0])):
            self.DataTable.setItem(
                row, 17, QtWidgets.QTableWidgetItem(alistoflists[0][row]))
            self.DataTable.setItem(
                row, 18, QtWidgets.QTableWidgetItem(alistoflists[1][row]))
            self.DataTable.setItem(
                row, 19, QtWidgets.QTableWidgetItem(alistoflists[2][row]))
            self.DataTable.setItem(
                row, 20, QtWidgets.QTableWidgetItem(alistoflists[3][row]))
            self.DataTable.setItem(
                row, 21, QtWidgets.QTableWidgetItem(alistoflists[4][row]))

    def Tablemaker3(self, alist):
        self.DataTable.setRowCount(len(alist))
        for row in range(len(alist)):
            self.DataTable.setItem(
                row, 22, QtWidgets.QTableWidgetItem(str(alist[row]['heterozygousALL'])))
            self.DataTable.setItem(
                row, 23, QtWidgets.QTableWidgetItem(str(alist[row]['minorhomozygousALL'])))
            self.DataTable.setItem(
                row, 24, QtWidgets.QTableWidgetItem(str(alist[row]['majorhomozygousALL'])))
            self.DataTable.setItem(
                row, 25, QtWidgets.QTableWidgetItem(str(alist[row]['heterozygousAFR'])))
            self.DataTable.setItem(
                row, 26, QtWidgets.QTableWidgetItem(str(alist[row]['minorhomozygousAFR'])))
            self.DataTable.setItem(
                row, 27, QtWidgets.QTableWidgetItem(str(alist[row]['majorhomozygousAFR'])))
            self.DataTable.setItem(
                row, 28, QtWidgets.QTableWidgetItem(str(alist[row]['heterozygousEUR'])))
            self.DataTable.setItem(
                row, 29, QtWidgets.QTableWidgetItem(str(alist[row]['minorhomozygousEUR'])))
            self.DataTable.setItem(
                row, 30, QtWidgets.QTableWidgetItem(str(alist[row]['majorhomozygousEUR'])))
            self.DataTable.setItem(
                row, 31, QtWidgets.QTableWidgetItem(str(alist[row]['heterozygousAMR'])))
            self.DataTable.setItem(
                row, 32, QtWidgets.QTableWidgetItem(str(alist[row]['minorhomozygousAMR'])))
            self.DataTable.setItem(
                row, 33, QtWidgets.QTableWidgetItem(str(alist[row]['majorhomozygousAMR'])))
            self.DataTable.setItem(
                row, 34, QtWidgets.QTableWidgetItem(str(alist[row]['heterozygousEAS'])))
            self.DataTable.setItem(
                row, 35, QtWidgets.QTableWidgetItem(str(alist[row]['minorhomozygousEAS'])))
            self.DataTable.setItem(
                row, 36, QtWidgets.QTableWidgetItem(str(alist[row]['majorhomozygousEAS'])))
            self.DataTable.setItem(
                row, 37, QtWidgets.QTableWidgetItem(str(alist[row]['heterozygousSAS'])))
            self.DataTable.setItem(
                row, 38, QtWidgets.QTableWidgetItem(str(alist[row]['minorhomozygousSAS'])))
            self.DataTable.setItem(
                row, 39, QtWidgets.QTableWidgetItem(str(alist[row]['majorhomozygousSAS'])))

        # Ui update while a thread starts/finishes
        #########################################################################

    def infosearch(self):
        Text = 'Fetching Summary\n\nPlease Wait...'
        self.label_8.setText(Text)

    def Ui_Load_UIDs_off(self):
        self.process.setText('Retrieving SNPs from DbSNP please wait...')
        self.Retmax.setEnabled(False)
        self.Retstart.setEnabled(False)
        self.UIDlist.setEnabled(False)
        self.RetrieveUIDs.setEnabled(False)
        self.Savebut.setEnabled(False)
        self.Clearbut.setEnabled(False)
        self.genename.setEnabled(False)
        self.clinsignificance.setEnabled(False)

    def Ui_Load_UIDs_on(self):
        self.process.setText('No process running: waiting for input')
        self.genename.setText('(optional)')
        self.Retmax.setEnabled(True)
        self.Retstart.setEnabled(True)
        self.UIDlist.setEnabled(True)
        self.commonfill.setEnabled(True)
        self.commonfill.setChecked(False)
        self.RetrieveUIDs.setEnabled(True)
        self.SelectAll.setEnabled(True)
        self.Savebut.setEnabled(True)
        self.Clearbut.setEnabled(True)
        self.genename.setEnabled(True)
        self.clinsignificance.setEnabled(True)
        self.Retmax.setValue(0)
        self.Retstart.setValue(0)

    def Ui_Get_Freqs_off(self):
        self.process.setText(
            'Retrieving SNP data please wait...')
        self.GetFreq.setEnabled(False)
        self.Savebut.setEnabled(False)
        self.Clearbut.setEnabled(False)
        self.DataTable.setEnabled(False)

    def Ui_Get_Freqs_on(self):
        self.process.setText('No process running: waiting for input')
        self.GetFreq.setEnabled(True)
        self.Savebut.setEnabled(True)
        self.Clearbut.setEnabled(True)
        self.DataTable.setEnabled(True)

    # Function section
    ###################################################################################################################

    # Function to save data from the table as a csv
    ########################################################################################################################################################

    def save_clicked(self):
        columnHeaders = []

        # create column header list
        for j in range(self.DataTable.model().columnCount()):
            columnHeaders.append(self.DataTable.horizontalHeaderItem(j).text())

        df = pd.DataFrame(columns=columnHeaders)

        # create dataframe object recordset
        for row in range(self.DataTable.rowCount()):
            for col in range(self.DataTable.columnCount()):
                df.at[row, columnHeaders[col]] = self.DataTable.item(
                    row, col).text()

        text, ok = QtWidgets.QInputDialog.getText(
            self, 'Saving csv', 'File Name:')

        if ok and text != '':
            textls = text.split('.')
            if len(textls) == 2:
                if textls[0].isalnum() == True and textls[1].isalnum() == True:
                    if textls[1] == 'csv':
                        df.to_csv(f'{text[0]}.{textls[1]}', index=False)
                        msg = QtWidgets.QMessageBox.information(self, ' ', 'File Saved.',
                                                                QtWidgets.QMessageBox.Ok)
                    else:
                        df.to_csv(f'{textls[0]}.csv', index=False)
                        msg = QtWidgets.QMessageBox.information(self, ' ', 'File Saved.',
                                                                QtWidgets.QMessageBox.Ok)
                else:
                    msg = QtWidgets.QMessageBox.warning(self, 'ERROR!', 'Only use alphanumeric characters\nfor the file name.',
                                                        QtWidgets.QMessageBox.Ok)

            elif len(textls) == 1 and text.isalnum() == True:
                df.to_csv(f'{text}.csv', index=False)
                msg = QtWidgets.QMessageBox.information(self, ' ', 'File Saved.',
                                                        QtWidgets.QMessageBox.Ok)
            else:
                msg = QtWidgets.QMessageBox.warning(self, 'ERROR!', 'Only use alphanumeric characters\nfor the file name.',
                                                    QtWidgets.QMessageBox.Ok)

        else:
            msg = QtWidgets.QMessageBox.warning(self, 'ERROR', 'A file name cant be blank',
                                                QtWidgets.QMessageBox.Ok)

 # Functions for searching and fetching SNP data from dbSNP and Ensemble
            ###############################################################################################################################################################################################

    def available_SNV(self, retstart, retmax, gene, clinsignificance):
        UIDlist = []
        '''Input = from wich row in the database should the results begin, the number of results you want
                Output = the list of available UIDs based on the inputs'''
        if clinsignificance == 'No Filtering':
            if gene == '(optional)' and self.commonfill.isChecked() == True:
                r = requests.get(
                    f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+AND+00000.0100:+00001.0000[GLOBAL_MAF]&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
            elif gene == '(optional)' and self.commonfill.isChecked() == False:
                r = requests.get(
                    f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+NOT+00000.0000[Global Minor Allele Frequency]&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
            elif gene != '(optional)' and self.commonfill.isChecked() == True:
                if gene.isnumeric() or (gene == 'Y' or gene == 'X'):
                    r = requests.get(
                        f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+AND+00000.0100:+00001.0000[GLOBAL_MAF]+AND+{gene}[Chromosome]&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
                elif (':' not in gene) and (not gene.isnumeric()) and (not(gene == 'Y' or gene == 'X')):
                    r = requests.get(
                        f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+AND+00000.0100:+00001.0000[GLOBAL_MAF]+AND+{gene}[Gene Name]&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
                elif (':' in gene) and ('-' in gene):
                    x = gene.split(':')
                    y = x[1].split('-')
                    chrom = x[0]
                    bpstart = y[0]
                    bpend = y[1]
                    r = requests.get(
                        f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+AND+00000.0100:+00001.0000[GLOBAL_MAF]+AND+({chrom}[Chromosome]+AND+({bpstart}[CHRPOS]+:+{bpend}[CHRPOS]))&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')

            elif gene != '(optional)' and self.commonfill.isChecked() == False:
                if gene.isnumeric() or (gene == 'Y' or gene == 'X'):
                    r = requests.get(
                        f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+NOT+00000.0000[Global Minor Allele Frequency]+AND+{gene}[Chromosome]&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
                elif (':' not in gene) and (not gene.isnumeric()) and (not(gene == 'Y' or gene == 'X')):
                    r = requests.get(
                        f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+NOT+00000.0000[Global Minor Allele Frequency]+AND+{gene}[Gene Name]&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
                elif (':' in gene) and ('-' in gene):
                    x = gene.split(':')
                    y = x[1].split('-')
                    chrom = x[0]
                    bpstart = y[0]
                    bpend = y[1]
                    r = requests.get(
                        f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+NOT+00000.0000[Global Minor Allele Frequency]+AND+{chrom}[Chromosome]+AND+{bpstart}[CHRPOS]+:+{bpend}[CHRPOS]&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
            data = json.loads(r.content)
            strlist = (data['esearchresult']['idlist'])
            numlist = [int(x) for x in strlist]
            UIDlist = sorted(numlist)
            self.label.setText(
                f"SNPs: {len(UIDlist)}/{data['esearchresult']['count']}")

        else:
            clinsignificance1 = str(clinsignificance).replace(' ', '+')
            if gene == '(optional)' and self.commonfill.isChecked() == True:
                r = requests.get(
                    f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+AND+00000.0100:+00001.0000[GLOBAL_MAF]+AND+{clinsignificance1}[Clinical Significance]&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
            elif gene == '(optional)' and self.commonfill.isChecked() == False:
                r = requests.get(
                    f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+NOT+00000.0000[Global Minor Allele Frequency]+AND+{clinsignificance1}[Clinical Significance]&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
            elif gene != '(optional)' and self.commonfill.isChecked() == True:
                if gene.isnumeric() or (gene == 'Y' or gene == 'X'):
                    r = requests.get(
                        f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+AND+00000.0100:+00001.0000[GLOBAL_MAF]+AND+{clinsignificance1}+AND+{gene}[Chromosome]&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
                elif (':' not in gene) and (not gene.isnumeric()) and (not(gene == 'Y' or gene == 'X')):
                    r = requests.get(
                        f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+AND+00000.0100:+00001.0000[GLOBAL_MAF]+AND+{clinsignificance1}+AND+{gene}[Gene Name]&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
                elif (':' in gene) and ('-' in gene):
                    x = gene.split(':')
                    y = x[1].split('-')
                    chrom = x[0]
                    bpstart = y[0]
                    bpend = y[1]
                    r = requests.get(
                        f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+AND+00000.0100:+00001.0000[GLOBAL_MAF]AND+{clinsignificance1}+AND+({chrom}[Chromosome]+AND+({bpstart}[CHRPOS]+:+{bpend}[CHRPOS]))&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
            elif gene != '(optional)' and self.commonfill.isChecked() == False:
                if gene.isnumeric() or (gene == 'Y' or gene == 'X'):
                    r = requests.get(
                        f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+NOT+00000.0000[Global Minor Allele Frequency]+AND+{clinsignificance1}[Clinical Significance]+AND+{gene}[Chromosome]&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
                elif (':' not in gene) and (not gene.isnumeric()) and (not(gene == 'Y' or gene == 'X')):
                    r = requests.get(
                        f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+NOT+00000.0000[Global Minor Allele Frequency]+AND+{clinsignificance1}[Clinical Significance]+AND+{gene}[Gene Name]&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
                elif (':' in gene) and ('-' in gene):
                    x = gene.split(':')
                    y = x[1].split('-')
                    chrom = x[0]
                    bpstart = y[0]
                    bpend = y[1]
                    r = requests.get(
                        f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=1000genomes+has+frequency+filter[Filter]+AND+snv[SNP Class]+NOT+00000.0000[Global Minor Allele Frequency]+AND+{clinsignificance1}[Clinical Significance]+AND+({chrom}[Chromosome]+AND+({bpstart}[CHRPOS]+:+{bpend}[CHRPOS]))&retstart={retstart}&retmax={retmax}&retmode=json&sort=SNP_ID')
            data = json.loads(r.content)
            strlist = (data['esearchresult']['idlist'])
            numlist = [int(x) for x in strlist]
            UIDlist = sorted(numlist)
            self.label.setText(
                f"SNPs: {len(UIDlist)}/{data['esearchresult']['count']}")

        return UIDlist

    def infosum(self, uid):
        '''creates the data summary for the SNV'''
        rsSNV = 'rs' + str(uid.text())
        server = "https://rest.ensembl.org"
        ext = f"/variation/human/{rsSNV}?phenotypes=1"
        r = requests.get(
            server+ext, headers={"Content-Type": "application/json"}).json()
        Seq = r['mappings'][0]['allele_string']
        Minor_Allele = r['minor_allele']
        func = r['most_severe_consequence']
        clinical = ''
        Traits = ''
        AltTraits = ''
        Gene = ''
        if r['phenotypes']:
            for items in r['phenotypes']:
                if 'risk_allele' in items:
                    if 'genes' in items:
                        Gene = items['genes']
                    if items['risk_allele'] == Minor_Allele:
                        Trait = items['trait']
                        if Trait not in Traits:
                            Traits += f'{Trait}\n'
                    elif items['risk_allele'] != Minor_Allele:
                        AltTrait = str(items['risk_allele']) + \
                            '->'+str(items['trait']).lower()
                        if AltTrait not in AltTraits:
                            AltTraits += f'{AltTrait}\n'
        try:
            for item in r['clinical_significance']:
                clinical += f'{item}\n'
        except KeyError:
            clinical += 'N\A'
        self.label_8.setText(f"Name={rsSNV}\n"
                             "\n"
                             f"Minor Allele={Minor_Allele}\n"
                             "\n"
                             f"Sequence={Seq}\n"
                             "\n"
                             f"Gene={Gene}\n"
                             "\n"
                             f"Function=\n{func}\n"
                             "\n"
                             f"Traits=\n{Traits.lower()}"
                             "\n"
                             f"Major Allele Traits=\n{AltTraits}"
                             "\n"
                             f"Clinical Significance=\n{clinical}"
                             )

    def Phenfinder(self, UIDlist):
        '''Finds the data to fill the table'''
        newls = [("rs" + str(n)) for n in UIDlist]
        UIDliststr = str(newls).replace("'", '"')
        server = "https://rest.ensembl.org"
        ext = "/variation/homo_sapiens?phenotypes=1"
        data = '{"ids":' + UIDliststr + '}'
        headers = {"Content-Type": "application/json",
                   "Accept": "application/json"}
        r = requests.post(server+ext, headers=headers,
                          data=data)

        decoded = json.loads(r.content)
        decodeduidsorted = sorted([int(uid.strip('rs'))
                                   for uid in decoded.keys()])
        true_uidlist = ['rs' + str(uid) for uid in decodeduidsorted]
        clinicallist = []
        funclist = []
        AltTraitlist = []
        Traitlist = []
        Genelist = []
        for uid in true_uidlist:
            AltTraitsstring = ''
            Traitstring = ''
            AltTraitsstring = ''
            Gene = ''
            try:
                funclist.append(decoded[uid]['most_severe_consequence'])
            except KeyError:
                funclist.append('Not specified')
            try:
                clinicallist.append(
                    str(decoded[uid]["clinical_significance"]).strip('[]'))
            except KeyError:
                clinicallist.append("Not specified")
            for items in decoded[uid]["phenotypes"]:
                try:
                    Gene = items['genes']
                except KeyError:
                    Gene = 'Not specified'
                try:
                    if items["risk_allele"] == decoded[uid]["minor_allele"]:
                        try:
                            Trait = items['trait']
                            if Trait.lower() not in Traitstring:
                                Traitstring += f'|{Trait.lower()}|'
                        except KeyError:
                            Traitstring = 'Not specified'
                    else:
                        try:
                            AltTrait = items['trait']
                            if AltTrait.lower() not in AltTraitsstring:
                                AltTraitsstring = f'|{AltTrait.lower()}|'
                        except KeyError:
                            AltTraitsstring = 'Not specified'
                except KeyError:
                    continue
            AltTraitlist.append(AltTraitsstring)
            Traitlist.append(Traitstring)
            Genelist.append(Gene)

        return [funclist, Genelist, Traitlist, AltTraitlist, clinicallist]

    def Popfinder(self, UIDlist):
        '''Finds the data to fill the table'''
        newls = [("rs" + str(n)) for n in UIDlist]
        UIDliststr = str(newls).replace("'", '"')
        server = "https://rest.ensembl.org"
        ext = "/variation/homo_sapiens?pops=1"
        data = '{"ids":' + UIDliststr + '}'
        headers = {"Content-Type": "application/json",
                   "Accept": "application/json"}
        r = requests.post(server+ext, headers=headers,
                          data=data)

        decoded = json.loads(r.content)
        chromosome = []
        popdata = []
        position = []
        minor_allele = []
        major_allele = []

        decodeduidsorted = sorted([int(uid.strip('rs'))
                                   for uid in decoded.keys()])
        true_uidlist = ['rs' + str(uid) for uid in decodeduidsorted]
        for uid in true_uidlist:
            popdict = {}
            Allslist = []
            AFRlist = []
            AMRlist = []
            EASlist = []
            EURlist = []
            SASlist = []
            major = ''
            minor = ''

            try:
                minor += decoded[uid]["minor_allele"]
                minor_allele.append(minor)
            except TypeError:
                minor_allele.append('NA')

            if len(decoded[uid]["mappings"]) > 1:
                try:
                    chromosome.append(
                        (decoded[uid]["mappings"][0]["location"].split(':'))[0])
                    position.append(
                        (decoded[uid]["mappings"][0]["location"].split('-'))[1])
                    major = decoded[uid]["mappings"][0]["allele_string"].split(
                        '/')[0]
                    if major != minor:
                        major_allele.append(major)
                        major = ''
                    else:
                        major_allele.append(decoded[uid]["mappings"][0]["allele_string"].split(
                            '/')[1])
                        major = ''

                except KeyError:
                    chromosome.append('NA')
                    position.append('NA')
            else:
                for items in decoded[uid]["mappings"]:
                    try:
                        chromosome.append((items["location"].split(':'))[0])
                        position.append(
                            (items["location"].split('-'))[1])
                        if len(items["allele_string"]) > 3:
                            major = items["ancestral_allele"]
                            major_allele.append(major)
                        else:
                            major = items["allele_string"].split('/')[0]
                            if major != minor:
                                major_allele.append(major)
                            else:
                                major_allele.append(
                                    items["allele_string"].split('/')[1])
                    except KeyError:
                        chromosome.append('NA')
                        position.append('NA')
            for items in decoded[uid]['populations']:
                if items['population'] == '1000GENOMES:phase_3:ALL':
                    Allslist.append(items["frequency"])
                    popdict['1000GENOMES:phase_3:ALL_major'] = max(Allslist)
                    if min(Allslist) == 1:
                        popdict['1000GENOMES:phase_3:ALL_minor'] = 0
                    else:
                        popdict['1000GENOMES:phase_3:ALL_minor'] = min(
                            Allslist)
                if items['population'] == '1000GENOMES:phase_3:AFR':
                    AFRlist.append(items["frequency"])
                    popdict['1000GENOMES:phase_3:AFR_major'] = max(AFRlist)
                    if min(AFRlist) == 1:
                        popdict['1000GENOMES:phase_3:AFR_minor'] = 0
                    else:
                        popdict['1000GENOMES:phase_3:AFR_minor'] = min(AFRlist)
                if items['population'] == '1000GENOMES:phase_3:AMR':
                    AMRlist.append(items["frequency"])
                    popdict['1000GENOMES:phase_3:AMR_major'] = max(AMRlist)
                    if min(AMRlist) == 1:
                        popdict['1000GENOMES:phase_3:AMR_minor'] = 0
                    else:
                        popdict['1000GENOMES:phase_3:AMR_minor'] = min(AMRlist)
                if items['population'] == '1000GENOMES:phase_3:EAS':
                    EASlist.append(items["frequency"])
                    popdict['1000GENOMES:phase_3:EAS_major'] = max(EASlist)
                    if min(EASlist) == 1:
                        popdict['1000GENOMES:phase_3:EAS_minor'] = 0
                    else:
                        popdict['1000GENOMES:phase_3:EAS_minor'] = min(EASlist)
                if items['population'] == '1000GENOMES:phase_3:EUR':
                    EURlist.append(items["frequency"])
                    popdict['1000GENOMES:phase_3:EUR_major'] = max(EURlist)
                    if min(EURlist) == 1:
                        popdict['1000GENOMES:phase_3:EUR_minor'] = 0
                    else:
                        popdict['1000GENOMES:phase_3:EUR_minor'] = min(EURlist)
                if items['population'] == '1000GENOMES:phase_3:SAS':
                    SASlist.append(items["frequency"])
                    popdict['1000GENOMES:phase_3:SAS_major'] = max(SASlist)
                    if min(SASlist) == 1:
                        popdict['1000GENOMES:phase_3:SAS_minor'] = 0
                    else:
                        popdict['1000GENOMES:phase_3:SAS_minor'] = min(SASlist)
            popdata.append(popdict)

        return [newls, chromosome, position, minor_allele, major_allele, popdata]

    def Genotypefinder(self, UIDlist):
        '''Finds the data to fill the table'''
        newls = [("rs" + str(n)) for n in UIDlist]
        UIDliststr = str(newls).replace("'", '"')
        server = "https://rest.ensembl.org"
        ext = "/variation/homo_sapiens?population_genotypes=1"
        data = '{"ids":' + UIDliststr + '}'
        headers = {"Content-Type": "application/json",
                   "Accept": "application/json"}
        r = requests.post(server+ext, headers=headers,
                          data=data)

        decoded = json.loads(r.content)
        decodeduidsorted = sorted([int(uid.strip('rs'))
                                   for uid in decoded.keys()])
        true_uidlist = ['rs' + str(uid) for uid in decodeduidsorted]

        keys = ['majorhomozygousALL', 'majorhomozygousAFR',
                'majorhomozygousAMR', 'majorhomozygousEAS',
                'majorhomozygousEUR', 'majorhomozygousSAS',
                'minorhomozygousALL', 'minorhomozygousAFR',
                'minorhomozygousAMR', 'minorhomozygousEAS',
                'minorhomozygousEUR', 'minorhomozygousSAS',
                'heterozygousALL', 'heterozygousAFR',
                'heterozygousAMR', 'heterozygousEAS',
                'heterozygousEUR', 'heterozygousSAS']

        Genotypelist = []

        for uid in true_uidlist:
            gendict = {}
            minorallele = decoded[uid]['minor_allele']

            for items in decoded[uid]['population_genotypes']:
                if items['population'] == '1000GENOMES:phase_3:ALL':
                    # filter heterozygotous genotypes
                    if (items['genotype'].split('|')[0] == minorallele or items['genotype'].split('|')[1] == minorallele) and not items['genotype'] == f'{minorallele}|{minorallele}':
                        gendict['heterozygousALL'] = items['frequency']
                    # filter Minor allele Homozygous genotypes
                    if items['genotype'] == f'{minorallele}|{minorallele}':
                        gendict['minorhomozygousALL'] = items['frequency']
                    # filter Major allele Homozygous genotypes
                    if (items['genotype'].split('|')[0] != minorallele) and (items['genotype'].split('|')[1] != minorallele):
                        gendict['majorhomozygousALL'] = items['frequency']
                if items['population'] == '1000GENOMES:phase_3:AFR':
                    if (items['genotype'].split('|')[0] == minorallele or items['genotype'].split('|')[1] == minorallele) and not items['genotype'] == f'{minorallele}|{minorallele}':
                        gendict['heterozygousAFR'] = items['frequency']
                    if items['genotype'] == f'{minorallele}|{minorallele}':
                        gendict['minorhomozygousAFR'] = items['frequency']
                    if (items['genotype'].split('|')[0] != minorallele) and (items['genotype'].split('|')[1] != minorallele):
                        gendict['majorhomozygousAFR'] = items['frequency']
                if items['population'] == '1000GENOMES:phase_3:AMR':
                    if (items['genotype'].split('|')[0] == minorallele or items['genotype'].split('|')[1] == minorallele) and not items['genotype'] == f'{minorallele}|{minorallele}':
                        gendict['heterozygousAMR'] = items['frequency']
                    if items['genotype'] == f'{minorallele}|{minorallele}':
                        gendict['minorhomozygousAMR'] = items['frequency']
                    if (items['genotype'].split('|')[0] != minorallele) and (items['genotype'].split('|')[1] != minorallele):
                        gendict['majorhomozygousAMR'] = items['frequency']
                if items['population'] == '1000GENOMES:phase_3:EAS':
                    if (items['genotype'].split('|')[0] == minorallele or items['genotype'].split('|')[1] == minorallele) and not items['genotype'] == f'{minorallele}|{minorallele}':
                        gendict['heterozygousEAS'] = items['frequency']
                    if items['genotype'] == f'{minorallele}|{minorallele}':
                        gendict['minorhomozygousEAS'] = items['frequency']
                    if (items['genotype'].split('|')[0] != minorallele) and (items['genotype'].split('|')[1] != minorallele):
                        gendict['majorhomozygousEAS'] = items['frequency']
                if items['population'] == '1000GENOMES:phase_3:EUR':
                    if (items['genotype'].split('|')[0] == minorallele or items['genotype'].split('|')[1] == minorallele) and not items['genotype'] == f'{minorallele}|{minorallele}':
                        gendict['heterozygousEUR'] = items['frequency']
                    if items['genotype'] == f'{minorallele}|{minorallele}':
                        gendict['minorhomozygousEUR'] = items['frequency']
                    if (items['genotype'].split('|')[0] != minorallele) and (items['genotype'].split('|')[1] != minorallele):
                        gendict['majorhomozygousEUR'] = items['frequency']
                if items['population'] == '1000GENOMES:phase_3:SAS':
                    if (items['genotype'].split('|')[0] == minorallele or items['genotype'].split('|')[1] == minorallele) and not items['genotype'] == f'{minorallele}|{minorallele}':
                        gendict['heterozygousSAS'] = items['frequency']
                    if items['genotype'] == f'{minorallele}|{minorallele}':
                        gendict['minorhomozygousSAS'] = items['frequency']
                    if (items['genotype'].split('|')[0] != minorallele) and (items['genotype'].split('|')[1] != minorallele):
                        gendict['majorhomozygousSAS'] = items['frequency']
                for items in keys:
                    if items not in gendict.keys():
                        gendict[items] = 0
            Genotypelist.append(gendict)

        return Genotypelist
    # Functions regarding the main thread/UI
    #########################################################################################################################################################################################################

    def Selected_UID_list(self):
        '''Iterates over the Selected UID list'''
        selectedlist = []
        for i in range(self.SelectedUID.count()):
            selectedlist.append(self.SelectedUID.item(i).text())
        return selectedlist

    def clear_all(self):
        '''Cleards all data'''
        self.UIDlist.clear()
        self.SelectedUID.clear()
        self.DataTable.clearContents()
        self.DataTable.model().removeRows(0, self.DataTable.rowCount())
        self.genename.setText('(optional)')
        self.genename.setEnabled(True)
        self.clinsignificance.setEnabled(True)
        self.label_2.setText("Selected SNPs:")
        self.label.setText("SNPs: shown/hits")
        self.label_8.setText("Name=None\n"
                             "\n"
                             "Minor Allele=None\n"
                             "\n"
                             "Sequence=None\n"
                             "\n"
                             "Gene=None\n"
                             "\n"
                             "Function=None\n"
                             "\n"
                             "Traits=None\n"
                             "\n"
                             "Major Allele Traits=None\n"
                             "\n"
                             "Clinical Significance=None\n"
                             )

    def clear_selection(self):
        '''Clears the selection'''
        self.SelectedUID.clear()
        self.label_2.setText("Selected SNPs:")

    def remove_item(self, item):
        '''Removes the item clicked from the Selected UIDs list'''
        self.SelectedUID.takeItem(self.SelectedUID.row(item))
        self.label_2.setText(str(len(self.Selected_UID_list())) + ' UIDs')

    def UID_list(self):
        '''Iterates over the UID list'''
        uidlist = []
        for i in range(self.UIDlist.count()):
            uidlist.append(self.UIDlist.item(i).text())
        return uidlist

    def get_all(self):
        for item in self.UID_list():
            if item not in self.Selected_UID_list():
                self.SelectedUID.addItem(item)
        self.SelectAll.setCheckState(False)
        self.label_2.setText(str(len(self.Selected_UID_list())) + ' SNPs')

    def get_one(self, item):
        '''Passes one UID into Selected UIDs list if it isnt already there'''
        if item.text() not in self.Selected_UID_list():
            self.SelectedUID.addItem(item.text())
        self.label_2.setText(str(len(self.Selected_UID_list())) + ' SNPs')
# Ui Elements
    ##############################################################################################################################################################################

    def setupUi(self, MainWindow):

        # connection to thread
        self.threadpool = QtCore.QThreadPool()

        # Mainwindow
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(900, 613)
        MainWindow.setMaximumSize(QtCore.QSize(900, 613))
        MainWindow.setSizeIncrement(QtCore.QSize(0, 0))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(255, 193, 131))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active,
                         QtGui.QPalette.ToolTipBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 193, 131))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive,
                         QtGui.QPalette.ToolTipBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 193, 131))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled,
                         QtGui.QPalette.ToolTipBase, brush)
        MainWindow.setPalette(palette)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(-10, 0, 911, 591))
        self.tabWidget.setMovable(True)
        self.tabWidget.setObjectName("tabWidget")
        self.FreqfinderTab = QtWidgets.QWidget()
        self.FreqfinderTab.setObjectName("FreqfinderTab")
        self.widget = QtWidgets.QWidget(self.FreqfinderTab)
        self.widget.setGeometry(QtCore.QRect(10, -10, 901, 601))
        self.widget.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.widget.setAutoFillBackground(False)
        self.widget.setStyleSheet("QWidget#widget{\n"
                                  "background-color:qlineargradient(spread:pad, x1:0, y1:0, x2:0.971591, y2:1, stop:0 rgba(170, 255, 255, 255), stop:1 rgba(255, 255, 255, 255));}")
        self.widget.setObjectName("widget")
        self.gridLayoutWidget = QtWidgets.QWidget(self.widget)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(6, 100, 281, 401))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.UIDlayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.UIDlayout.setContentsMargins(0, 0, 0, 0)
        self.UIDlayout.setSpacing(10)
        self.UIDlayout.setObjectName("UIDlayout")
        self.label_4 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_4.setObjectName("label_4")
        self.UIDlayout.addWidget(self.label_4, 1, 1, 1, 1)

        # UIDlist
        self.UIDlist = QtWidgets.QListWidget(self.gridLayoutWidget)
        self.UIDlist.itemDoubleClicked.connect(self.get_one)
        self.UIDlist.itemClicked.connect(self.UIDclicked)
        self.UIDlist.setObjectName("UIDlist")
        self.UIDlayout.addWidget(self.UIDlist, 7, 1, 1, 2)

        # static
        self.label_10 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_10.setObjectName("label_10")
        self.UIDlayout.addWidget(self.label_10, 6, 4, 1, 2)
        self.line = QtWidgets.QFrame(self.gridLayoutWidget)
        self.line.setFrameShadow(QtWidgets.QFrame.Plain)
        self.line.setMidLineWidth(1)
        self.line.setFrameShape(QtWidgets.QFrame.VLine)
        self.line.setObjectName("line")
        self.UIDlayout.addWidget(self.line, 1, 3, 2, 1)
        self.label_3 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_3.setObjectName("label_3")
        self.UIDlayout.addWidget(self.label_3, 2, 1, 1, 1)

        # Retstart Spinbox
        self.Retstart = QtWidgets.QSpinBox(self.gridLayoutWidget)
        self.Retstart.setStatusTip("")
        self.Retstart.setRange(0, 38216813)
        self.Retstart.setObjectName("Retstart")
        self.UIDlayout.addWidget(self.Retstart, 1, 2, 1, 1)

        # Retrieve Uids Button
        self.RetrieveUIDs = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.RetrieveUIDs.clicked.connect(self.RetrieveSNPs_clicked)
        self.RetrieveUIDs.setObjectName("RetrieveUIDs")
        self.UIDlayout.addWidget(self.RetrieveUIDs, 5, 1, 1, 5)

        self.label = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label.setObjectName("label")
        self.UIDlayout.addWidget(self.label, 6, 1, 1, 2)

        # Retmax Spinbox
        self.Retmax = QtWidgets.QSpinBox(self.gridLayoutWidget)
        self.Retmax.setRange(0, 100)
        self.Retmax.setObjectName("Retmax")
        self.UIDlayout.addWidget(self.Retmax, 2, 2, 1, 1)

        self.commonfill = QtWidgets.QCheckBox(self.gridLayoutWidget)
        self.commonfill.setAcceptDrops(False)
        self.commonfill.setObjectName("commonfill")
        self.UIDlayout.addWidget(self.commonfill, 1, 4, 2, 2)

        self.label_5 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_5.setObjectName("label_5")
        self.UIDlayout.addWidget(self.label_5, 4, 1, 1, 2)

        # Gene name line edit
        self.genename = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.genename.setText('(optional)')
        self.genename.setReadOnly(False)
        self.genename.setObjectName("genename")
        self.UIDlayout.addWidget(self.genename, 3, 3, 1, 3)

        # clinsignificance line edit
        self.clinsignificance = QtWidgets.QComboBox(self.gridLayoutWidget)
        self.clinsignificance.addItems(['No Filtering', 'benign', 'benign likely benign', 'likely benign', 'pathogenic', 'pathogenic likely pathogenic', 'likely pathogenic',
                                        'conflicting interpretations of pathogenicity', 'drug response', 'protective', 'risk factor', 'uncertain significance'])
        self.clinsignificance.setObjectName("clinsignificance")

        self.UIDlayout.addWidget(self.clinsignificance, 4, 3, 1, 3)

        self.label_11 = QtWidgets.QLabel(self.gridLayoutWidget)
        font = QtGui.QFont()
        font.setPointSize(8)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.UIDlayout.addWidget(self.label_11, 3, 1, 1, 2)

        # Summary Label/scroll area
        self.scrollArea = QtWidgets.QScrollArea(self.gridLayoutWidget)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 125, 225))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(
            self.scrollAreaWidgetContents)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_8 = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.label_8.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.label_8.setObjectName("label_8")
        self.verticalLayout_2.addWidget(self.label_8)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.UIDlayout.addWidget(self.scrollArea, 7, 3, 1, 3)

        # Static
        self.gridLayoutWidget_2 = QtWidgets.QWidget(self.widget)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(420, 100, 471, 401))
        self.gridLayoutWidget_2.setObjectName("gridLayoutWidget_2")
        self.Freqlayout = QtWidgets.QGridLayout(self.gridLayoutWidget_2)
        self.Freqlayout.setContentsMargins(0, 0, 0, 0)
        self.Freqlayout.setObjectName("Freqlayout")
        self.label_7 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_7.setStyleSheet("text-align: center;")
        self.label_7.setObjectName("label_7")
        self.Freqlayout.addWidget(self.label_7, 0, 0, 1, 1)

        # Freq Table
        self.DataTable = QtWidgets.QTableWidget(self.gridLayoutWidget_2)
        self.DataTable.setEnabled(True)
        self.DataTable.setColumnCount(40)
        self.DataTable.setHorizontalHeaderLabels(
            [' SNP ', 'Chromosome', 'Position', 'Minor allele', 'Major allele', 'Total minor allele frequency', 'Total major allele frequency', 'African minor allele frequency ',
             'African major allele frequency', 'European minor allele frequency', 'European major allele frequency', 'American minor allele frequency', 'American major allele frequency',
             'East Asian minor allele frequency', 'East Asian major allele frequency', 'South Asian minor allele frequency', 'South Asian major allele frequency',
             'Function', 'Gene', 'Minor allele traits', 'Major allele traits', 'Clinical Significance', 'Total heterozygous', 'Total minor allele homozygous', 'Total major allele homozygous',
             'African heterozygous', 'African minor allele homozygous', 'African major allele homozygous', 'European heterozygous', 'European minor allele homozygous', 'European major allele homozygous',
             'American heterozygous', 'American minor allele homozygous', 'American major allele homozygous', 'East Asian heterozygous', 'East Asian minor allele homozygous', 'East Asian major allele homozygous',
             'South Asian heterozygous', 'South Asian minor allele homozygous', 'South Asian major allele homozygous'])
        self.DataTable.resizeColumnsToContents()
        self.DataTable.resizeRowsToContents()
        self.DataTable.setObjectName("DataTable")
        self.Freqlayout.addWidget(self.DataTable, 6, 0, 1, 1)

        self.newhorizlayout = QtWidgets.QHBoxLayout()
        self.newhorizlayout.setObjectName("newhorizlayout")

        # GetfreqButton
        self.GetFreq = QtWidgets.QPushButton(self.gridLayoutWidget_2)
        self.GetFreq.setObjectName("GetFreq")
        self.GetFreq.clicked.connect(self.Get_Data_Clicked)
        self.newhorizlayout.addWidget(self.GetFreq)
        self.Freqlayout.addLayout(self.newhorizlayout, 2, 0, 1, 1)

        # static
        self.line_5 = QtWidgets.QFrame(self.gridLayoutWidget_2)
        self.line_5.setFrameShadow(QtWidgets.QFrame.Plain)
        self.line_5.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_5.setObjectName("line_5")
        self.Freqlayout.addWidget(self.line_5, 1, 0, 1, 1)
        self.verticalLayoutWidget = QtWidgets.QWidget(self.widget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(0, 0, 870, 110))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setSpacing(0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.Title = QtWidgets.QLabel(self.verticalLayoutWidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(24)
        font.setBold(False)
        font.setItalic(True)
        font.setUnderline(False)
        font.setWeight(50)
        self.Title.setFont(font)
        self.Title.setStyleSheet("")
        self.Title.setObjectName("Title")
        self.verticalLayout.addWidget(self.Title)
        self.line_2 = QtWidgets.QFrame(self.verticalLayoutWidget)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Plain)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setObjectName("line_2")
        self.verticalLayout.addWidget(self.line_2)
        self.label_6 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_6.setObjectName("label_6")
        self.verticalLayout.addWidget(self.label_6)
        self.gridLayoutWidget_3 = QtWidgets.QWidget(self.widget)
        self.gridLayoutWidget_3.setGeometry(QtCore.QRect(10, 500, 881, 73))
        self.gridLayoutWidget_3.setObjectName("gridLayoutWidget_3")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget_3)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        # Save Button
        self.Savebut = QtWidgets.QPushButton(self.gridLayoutWidget_3)
        self.Savebut.setSizeIncrement(QtCore.QSize(0, 0))
        self.Savebut.clicked.connect(self.save_clicked)
        self.Savebut.setObjectName("Savebut")
        self.gridLayout.addWidget(self.Savebut, 0, 3, 1, 1)

        # Clear all button
        self.Clearbut = QtWidgets.QPushButton(self.gridLayoutWidget_3)
        self.Clearbut.clicked.connect(self.clear_all)
        self.Clearbut.setObjectName("Clearbut")
        self.gridLayout.addWidget(self.Clearbut, 0, 2, 1, 1)

        # Static
        self.line_3 = QtWidgets.QFrame(self.gridLayoutWidget_3)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Plain)
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setObjectName("line_3")
        self.gridLayout.addWidget(self.line_3, 0, 1, 1, 1)

        # process label
        self.process = QtWidgets.QLabel(self.gridLayoutWidget_3)
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(True)
        font.setUnderline(False)
        font.setWeight(75)
        self.process.setFont(font)
        self.process.setObjectName("process")
        self.gridLayout.addWidget(self.process, 0, 0, 1, 1)

        # static
        self.gridLayoutWidget_4 = QtWidgets.QWidget(self.widget)
        self.gridLayoutWidget_4.setGeometry(QtCore.QRect(300, 100, 91, 401))
        self.gridLayoutWidget_4.setObjectName("gridLayoutWidget_4")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.gridLayoutWidget_4)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setHorizontalSpacing(10)
        self.gridLayout_2.setVerticalSpacing(12)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setSpacing(6)
        self.horizontalLayout.setObjectName("horizontalLayout")

        # Select All button
        self.SelectAll = QtWidgets.QCheckBox(self.gridLayoutWidget_4)
        self.SelectAll.stateChanged.connect(self.get_all)
        self.SelectAll.setObjectName("SelectAll")
        self.horizontalLayout.addWidget(self.SelectAll)

        # static
        self.gridLayout_2.addLayout(self.horizontalLayout, 0, 1, 1, 1)

        # selected label
        self.label_2 = QtWidgets.QLabel(self.gridLayoutWidget_4)
        self.label_2.setObjectName("label_2")

        self.gridLayout_2.addWidget(self.label_2, 2, 1, 1, 1)
        self.line_4 = QtWidgets.QFrame(self.gridLayoutWidget_4)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Plain)
        self.line_4.setMidLineWidth(1)
        self.line_4.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_4.setObjectName("line_4")
        self.gridLayout_2.addWidget(self.line_4, 0, 0, 4, 1)

        # clear selection button
        self.Clearselection = QtWidgets.QPushButton(self.gridLayoutWidget_4)
        self.Clearselection.clicked.connect(self.clear_selection)
        self.Clearselection.setObjectName("Clearselection")
        self.gridLayout_2.addWidget(self.Clearselection, 1, 1, 1, 1)

        # Selected UID list
        self.SelectedUID = QtWidgets.QListWidget(self.gridLayoutWidget_4)
        self.SelectedUID.setEnabled(True)
        self.SelectedUID.setMouseTracking(False)
        self.SelectedUID.itemClicked.connect(self.remove_item)
        self.SelectedUID.setObjectName("SelectedUID")
        self.gridLayout_2.addWidget(self.SelectedUID, 3, 1, 1, 1)

        # Second Tab
        self.tabWidget.addTab(self.FreqfinderTab, "")

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 900, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label_4.setText(_translate("MainWindow", "Start:"))
        self.label_5.setText(_translate(
            "MainWindow", "Clinical Significance\nSNP Filter:"))
        self.label_11.setText(_translate(
            "MainWindow", "Search only in a Gene:"))
        self.UIDlist.setToolTip(_translate(
            "MainWindow", "<html><head/><body><p align=\"center\">Click to view the SNP\'s info Double Click to add it to your selection<br/></p></body></html>"))
        self.label_10.setText(_translate(
            "MainWindow", " Summary:                     "))
        self.label_3.setText(_translate("MainWindow", "Stop:"))
        self.Retstart.setToolTip(_translate(
            "MainWindow", "<html><head/><body><p align=\"center\">Sequential index of the first SNP Retrieved</p></body></html>"))
        self.RetrieveUIDs.setToolTip(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
                                                "<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
                                                "p, li { white-space: pre-wrap; }\n"
                                                "</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:400; font-style:normal;\">\n"
                                                "<p align=\"center\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Retrieve Available SNPs from the DbSNP database </p>\n"
                                                "<p align=\"center\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">based on Filters and the Retstart and Retmax values.</p></body></html>"))
        self.RetrieveUIDs.setText(_translate("MainWindow", "Retrieve SNPs"))
        self.label.setText(_translate(
            "MainWindow", "SNPs: Shown/Hits"))
        self.Retmax.setToolTip(_translate(
            "MainWindow", "<html><head/><body><p align=\"center\">The number of SNPs to be Retrieved</p></body></html>"))
        self.genename.setToolTip(_translate(
            "MainWindow", "<html><head/><body><p>Chromosome: Enter a chromosome number number (e.g. '8')</p><p>Gene: Enter a gene ID (e.g. 'TP53')</p><p>Region: Enter a chromosome number and a genomic region (e.g. '8:1-50000') </p><</body></html>"))
        self.genename.setPlaceholderText(
            _translate("MainWindow", "e.g. 1"))
        self.label_11.setText(_translate(
            "MainWindow", "Search only in a Gene/\nChromosome/Region:"))
        self.label_8.setText(_translate("MainWindow", "Name=None\n"
                                        "\n"
                                        "Minor Allele=None\n"
                                        "\n"
                                        "Sequence=None\n"
                                        "\n"
                                        "Gene=None\n"
                                        "\n"
                                        "Function=None\n"
                                        "\n"
                                        "Traits=None\n"
                                        "\n"
                                        "Major Allele Traits=None\n"
                                        "\n"
                                        "Clinical Significance=None\n"
                                        ))
        self.commonfill.setToolTip(_translate(
            "MainWindow", "fill it to search for variants with a frequency range of (0.0100, 1.0000)"))
        self.commonfill.setText(_translate("MainWindow", " Search only\n"
                                           " for common\n"
                                           " variants "))
        self.label_7.setText(_translate(
            "MainWindow", "                                                       For the selected SNPs:"))
        self.GetFreq.setToolTip(_translate(
            "MainWindow", "<html><head/><body><p>Gets the frequency and other data for up to 100 SNPs</p></body></html>"))
        self.GetFreq.setText(_translate("MainWindow", "Get Data"))
        self.Title.setText(_translate("MainWindow", "SNP Finder"))
        self.label_6.setText(_translate(
            "MainWindow", "<html><head/><body><p><span style=\" font-size:11pt;\">Search databases for single nucleotide variations &amp; retrieve various data in a csv form.</span></p><p><br/></p></body></html>"))
        self.Savebut.setToolTip(_translate(
            "MainWindow", "<html><head/><body><p align=\"center\">Saves the data in csv file form </p></body></html>"))
        self.Savebut.setText(_translate("MainWindow", "Save data"))
        self.Clearbut.setToolTip(_translate(
            "MainWindow", "Clears all the data from the lists and the table"))
        self.Clearbut.setText(_translate(
            "MainWindow", "Clear all the data/Reset"))
        self.process.setText(_translate(
            "MainWindow", "No process running: waiting for input"))
        self.SelectAll.setToolTip(_translate(
            "MainWindow", "<html><head/><body><p align=\"center\">Selects all the SNPs Retrieved from the Available SNPs list</p></body></html>"))
        self.SelectAll.setText(_translate("MainWindow", "Select All"))
        self.label_2.setText(_translate("MainWindow", "Selected SNPs:"))
        self.Clearselection.setText(_translate("MainWindow", "Clear "))
        self.tabWidget.setTabText(self.tabWidget.indexOf(
            self.FreqfinderTab), _translate("MainWindow", "Allele Frequency Finder"))


# Application main
##########################
if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle('fusion')
    app.setWindowIcon(QtGui.QIcon('Appicon.png'))
    app.setApplicationName('Allele Frequency Finder')
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
