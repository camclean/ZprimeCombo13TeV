# -*- coding: utf-8 -*-
import datetime, os.path, shutil, config

class table:
   def __init__(self):
       self.columns = []
       self.pretty_colnames = []
       self.rows = []
       self.current_row = {}
       
   def add_column(self, colname, pretty_colname = None):
      if pretty_colname is None: pretty_colname = colname
      self.columns.append(colname)
      self.pretty_colnames.append(pretty_colname)
   
   # the html formatted stuff here:
   def set_column(self, colname, value):
       self.current_row[colname] = {'html': value}

   def get_columns(self):
       return self.columns
       
   def set_column_multiformat(self, colname, raw, **values):
       self.current_row[colname] = values
       self.current_row[colname]['raw'] = raw
   
   def add_row(self):
       row = []
       for c in self.columns: row.append(self.current_row[c])
       self.rows.append(row)
       self.current_row = {}
       
   def tex(self):
       result = r'\begin{tabular}{|'
       for c in self.pretty_colnames: result += "l|"
       result += "}\hline\n"
       result += " & ".join(self.pretty_colnames)
       result += r"\\" + "\n \hline \n"
       for r in self.rows:
           try:  result += ' & '.join([table._cellval(col, 'tex') for col in r])
           except KeyError: raise RuntimeError, "no tex content for column found"
           result += r"\\" + "\n"
       result += r"\hline" + "\n" +  r"\end{tabular}"
       return result

   def get_raw_rows(self):
       result = []
       for r in self.rows: result.append([cell['raw'] for cell in r])
       return result
       
   @staticmethod
   def _cellval(cell, format):
       if format in cell: return cell[format]
       else: return cell['raw']
       

   def html(self):
       result = '<table cellpadding="2" border="1">\n<tr>'
       for c in self.pretty_colnames:
           result += '<th>%s</th>' % c
       result += '</tr>\n'
       for r in self.rows:
          result += '<tr><td>'
          result += '</td><td>'.join([table._cellval(col, 'html') for col in r])
          result += '</td></tr>\n'
       result += '</table>\n'
       return result

class html_report:
    def __init__(self, output_filename):
        self.filename = output_filename
        self.f = open(output_filename, 'w')
        self.section_counter = 1
        self.section_active = False
        self.f.write(
        """<?xml version="1.0" ?>
    <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
    <html xmlns="http://www.w3.org/1999/xhtml">
    <head><title>theta-auto Report</title>
    <script type="text/javascript" src="js/prototype.js"></script>
    <script type="text/javascript" src="js/scriptaculous.js?load=effects,builder"></script>
    <script type="text/javascript" src="js/lightbox.js"></script>
    <link rel="stylesheet" href="css/lightbox.css" type="text/css" media="screen" />
    <style type="text/css">
    body {margin: 30px; font-family: Verdana, Arial, sans-serif; font-size: 12px; background:#eee;}
    h1, h2, h3 { padding: 0.2ex; padding-left:.8ex;}
    h1 {font-size: 150%%; background: #d44;}
    h2 {font-size: 110%%; background: #f88; padding-left:1.3ex;}
    h3 {font-size: 100%%; background: #aaa; padding-left:1.3ex;}
    .inner {padding-left: 3ex;}
    a img{ border: none;}
    p {margin:0 0 0 0; margin-top: 1.3ex;}
    ul {margin:0 0 0 0;}
    </style>
    </head><body>
    <p>Hint: click on top-level headers to toggle visibility of that section.</p>
    """)

    def new_section(self, heading, content = ''):
        if self.section_active: self.f.write('\n</div>\n')
        self.f.write('<h1 onclick=$(\'div%d\').toggle()>%d. %s</h1><div class="inner" id="div%d">\n' % (self.section_counter,
            self.section_counter, heading, self.section_counter))
        self.section_counter += 1
        self.f.write(content)
        self.section_active = True
        
    def add_p(self, text):
        self.f.write('<p>' + text + '</p>\n')
        
    def add_html(self, raw_html):
        self.f.write('<p>' + raw_html + '</p>\n')
        
    def write_html(self, html_target_path = ''):
        if self.section_active: self.f.write('\n</div>\n')
        global workdir, thetautils_dir
        d = {}
        d['datetime'] = datetime.datetime.now().isoformat(' ')
        d['workdir'] = config.workdir
        self.f.write("<hr /><p>This page was generated at %(datetime)s for workdir '%(workdir)s'.</p>" % d)
        self.f.write('</body></html>')
        self.f.close()
        if html_target_path!='':
            thetautils_dir = os.path.join(config.theta_dir, 'utils')
            if not os.path.exists(html_target_path):
                os.mkdir(html_target_path)
            if not os.path.isdir(html_target_path): raise RuntimeError, 'Given html_target_path "%s" is not a directory!' % html_target_path
            shutil.copy(self.filename, html_target_path)
            shutil.rmtree(os.path.join(html_target_path, 'js'), ignore_errors = True)
            shutil.rmtree(os.path.join(html_target_path, 'css'), ignore_errors = True)
            shutil.rmtree(os.path.join(html_target_path, 'images'), ignore_errors = True)
            shutil.copytree(os.path.join(thetautils_dir, 'htmlstuff', 'js'), os.path.join(html_target_path, 'js'))
            shutil.copytree(os.path.join(thetautils_dir, 'htmlstuff', 'css'), os.path.join(html_target_path, 'css'))
            shutil.copytree(os.path.join(thetautils_dir, 'htmlstuff', 'images'), os.path.join(html_target_path, 'images'))
            if os.path.exists(os.path.join(config.workdir, 'plots')):
                shutil.rmtree(os.path.join(html_target_path, 'plots'), ignore_errors = True)
                shutil.copytree(os.path.join(config.workdir, 'plots'), os.path.join(html_target_path, 'plots'))
        
